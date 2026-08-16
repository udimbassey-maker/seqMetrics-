import streamlit as st
from collections import Counter

from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils import MeltingTemp as mt


# Set up the web page layout, title, and browser icon
st.set_page_config(
    page_title="seqMetrics",
    page_icon="🧬",
    layout="centered"
)


# Standard genetic code dictionary mapping codons to their respective amino acids and properties
GENETIC_CODE = {
    "TTT": ("Phe", "Phenylalanine", ""),
    "TTC": ("Phe", "Phenylalanine", ""),
    "TTA": ("Leu", "Leucine", ""),
    "TTG": ("Leu", "Leucine", ""),

    "TCT": ("Ser", "Serine", ""),
    "TCC": ("Ser", "Serine", ""),
    "TCA": ("Ser", "Serine", ""),
    "TCG": ("Ser", "Serine", ""),

    "TAT": ("Tyr", "Tyrosine", ""),
    "TAC": ("Tyr", "Tyrosine", ""),
    "TAA": ("STOP", "Stop codon", "Stop codon"),
    "TAG": ("STOP", "Stop codon", "Stop codon"),

    "TGT": ("Cys", "Cysteine", ""),
    "TGC": ("Cys", "Cysteine", ""),
    "TGA": ("STOP", "Stop codon", "Stop codon"),
    "TGG": ("Trp", "Tryptophan", ""),

    "CTT": ("Leu", "Leucine", ""),
    "CTC": ("Leu", "Leucine", ""),
    "CTA": ("Leu", "Leucine", ""),
    "CTG": ("Leu", "Leucine", ""),

    "CCT": ("Pro", "Proline", ""),
    "CCC": ("Pro", "Proline", ""),
    "CCA": ("Pro", "Proline", ""),
    "CCG": ("Pro", "Proline", ""),

    "CAT": ("His", "Histidine", ""),
    "CAC": ("His", "Histidine", ""),
    "CAA": ("Gln", "Glutamine", ""),
    "CAG": ("Gln", "Glutamine", ""),

    "CGT": ("Arg", "Arginine", ""),
    "CGC": ("Arg", "Arginine", ""),
    "CGA": ("Arg", "Arginine", ""),
    "CGG": ("Arg", "Arginine", ""),

    "ATT": ("Ile", "Isoleucine", ""),
    "ATC": ("Ile", "Isoleucine", ""),
    "ATA": ("Ile", "Isoleucine", ""),
    "ATG": ("Met", "Methionine", "Start codon"),

    "ACT": ("Thr", "Threonine", ""),
    "ACC": ("Thr", "Threonine", ""),
    "ACA": ("Thr", "Threonine", ""),
    "ACG": ("Thr", "Threonine", ""),

    "AAT": ("Asn", "Asparagine", ""),
    "AAC": ("Asn", "Asparagine", ""),
    "AAA": ("Lys", "Lysine", ""),
    "AAG": ("Lys", "Lysine", ""),

    "AGT": ("Ser", "Serine", ""),
    "AGC": ("Ser", "Serine", ""),
    "AGA": ("Arg", "Arginine", ""),
    "AGG": ("Arg", "Arginine", ""),

    "GTT": ("Val", "Valine", ""),
    "GTC": ("Val", "Valine", ""),
    "GTA": ("Val", "Valine", ""),
    "GTG": ("Val", "Valine", ""),

    "GCT": ("Ala", "Alanine", ""),
    "GCC": ("Ala", "Alanine", ""),
    "GCA": ("Ala", "Alanine", ""),
    "GCG": ("Ala", "Alanine", ""),

    "GAT": ("Asp", "Aspartic acid", ""),
    "GAC": ("Asp", "Aspartic acid", ""),
    "GAA": ("Glu", "Glutamic acid", ""),
    "GAG": ("Glu", "Glutamic acid", ""),

    "GGT": ("Gly", "Glycine", ""),
    "GGC": ("Gly", "Glycine", ""),
    "GGA": ("Glycine", "Glycine", ""),
    "GGG": ("Gly", "Glycine", "")
}


# Remove any accidental whitespace and convert letters to uppercase
def clean_sequence(sequence):
    return "".join(sequence.upper().split())


# Check if the sequence contains only accepted standard nucleotides or ambiguous bases
def validate_sequence(sequence):
    valid_bases = set("ATGCN")
    return all(base in valid_bases for base in sequence)


# Calculate the percentage of Guanine and Cytosine in the sequence
def calculate_gc(sequence):
    if not sequence:
        return 0.0

    return gc_fraction(sequence, ambiguous="ignore") * 100


# Generate the reverse complement strand oriented 5' to 3'
def get_reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }

    return "".join(
        complement[base]
        for base in reversed(sequence)
    )


# Convert DNA coding strand into its corresponding RNA transcript by substituting Thymine with Uracil
def transcribe_coding_dna(sequence):
    return sequence.replace("T", "U")


# Count the exact occurrence of each nucleotide type in the sequence
def nucleotide_composition(sequence):
    counts = Counter(sequence)

    return {
        "A": counts.get("A", 0),
        "T": counts.get("T", 0),
        "G": counts.get("G", 0),
        "C": counts.get("C", 0),
        "N": counts.get("N", 0)
    }


# Split the sequence into consecutive 3-base triplets based on the selected reading frame
def get_codons(sequence, frame):
    return [
        sequence[i:i + 3]
        for i in range(frame, len(sequence) - 2, 3)
    ]


# Translate nucleotide triplets into an amino acid sequence
def translate_sequence(sequence, frame):
    codons = get_codons(sequence, frame)

    codon_details = []
    protein_sequence = []

    for codon in codons:

        if "N" in codon:
            amino_acid = "X"
            full_name = "Unknown amino acid"
            annotation = "Ambiguous codon"

        else:
            amino_acid, full_name, annotation = GENETIC_CODE[codon]

        codon_details.append({
            "Codon": codon,
            "Amino acid": amino_acid,
            "Name": full_name,
            "Annotation": annotation
        })

        if amino_acid == "STOP":
            protein_sequence.append("*")
            break

        protein_sequence.append(amino_acid)

    return codon_details, "".join(protein_sequence)


# Estimate the melting temperature based on DNA and salt concentrations
def calculate_melting_temperature(
    sequence,
    dna_concentration,
    salt_concentration
):
    if not sequence or "N" in sequence:
        return None

    return mt.Tm_NN(
        sequence,
        dnac=dna_concentration,
        Na=salt_concentration
    )


# Compute the total molecular mass of the single-stranded DNA sequence in Daltons
def calculate_molecular_mass(sequence):
    if not sequence or "N" in sequence:
        return None

    return molecular_weight(
        sequence,
        seq_type="DNA",
        double_stranded=False,
        circular=False,
        monoisotopic=False
    )


# Main layout and user interface design for the app
st.title("seqMetrics")

st.write(
    "DNA sequence analysis for sequence composition, "
    "GC content, transcription, reverse complement, "
    "translation, melting temperature, and molecular mass."
)


# Text box widget allowing users to input their custom DNA sequence
raw_sequence = st.text_area(
    "Enter DNA sequence (5' to 3')",
    placeholder="Example: ATGTACTGG"
)


if raw_sequence:

    sequence = clean_sequence(raw_sequence)

    if not validate_sequence(sequence):

        st.error(
            "Invalid DNA sequence. Use only A, T, G, C, or N."
        )

    else:

        st.success("DNA sequence validated successfully.")

        # Organize analysis metrics into clean, interactive application tabs
        tab1, tab2, tab3, tab4 = st.tabs(
            [
                "📊 Basic Metrics",
                "🔄 Sequence Operations",
                "🧬 Translation",
                "🧪 Physicochemical Properties"
            ]
        )


        # Tab 1: Display length, GC content, and base frequency breakdown
        with tab1:

            st.subheader("Basic Metrics")

            length = len(sequence)
            gc_content = calculate_gc(sequence)

            col1, col2 = st.columns(2)

            col1.metric(
                "Sequence Length",
                f"{length} bp"
            )

            col2.metric(
                "GC Content",
                f"{gc_content:.2f}%"
            )

            st.write("### Nucleotide Composition")

            composition = nucleotide_composition(sequence)

            col1, col2, col3, col4, col5 = st.columns(5)

            col1.metric("Adenine (A)", composition["A"])
            col2.metric("Thymine (T)", composition["T"])
            col3.metric("Guanine (G)", composition["G"])
            col4.metric("Cytosine (C)", composition["C"])
            col5.metric("Unknown (N)", composition["N"])

            if composition["N"] > 0:
                st.info(
                    "N represents an unknown or ambiguous nucleotide."
                )


        # Tab 2: Handle DNA-to-RNA transcription and reverse complement calculations
        with tab2:

            st.subheader("Sequence Operations")

            st.write("### DNA Sequence")

            st.code(sequence)

            st.write("### RNA Transcript")

            rna_sequence = transcribe_coding_dna(sequence)

            st.code(rna_sequence)

            st.caption(
                "The input is treated as a coding DNA strand. "
                "Transcription replaces T with U."
            )

            st.write("### Reverse Complement")

            reverse_complement = get_reverse_complement(sequence)

            st.code(reverse_complement)

            st.caption(
                "The reverse complement is reported 5' to 3'."
            )


        # Tab 3: Perform reading frame selection and protein sequence translation
        with tab3:

            st.subheader("DNA Translation")

            st.write(
                "The input sequence is treated as a coding DNA strand."
            )

            frame_option = st.selectbox(
                "Reading frame",
                [
                    "Frame +1",
                    "Frame +2",
                    "Frame +3"
                ]
            )

            frame_map = {
                "Frame +1": 0,
                "Frame +2": 1,
                "Frame +3": 2
            }

            frame = frame_map[frame_option]

            codon_details, protein_sequence = translate_sequence(
                sequence,
                frame
            )

            if not codon_details:

                st.warning(
                    "There are not enough nucleotides to form "
                    "a complete codon in this reading frame."
                )

            else:

                st.write("### Codon Translation")

                for detail in codon_details:

                    codon = detail["Codon"]
                    amino_acid = detail["Amino acid"]
                    full_name = detail["Name"]
                    annotation = detail["Annotation"]

                    if annotation == "Start codon":

                        st.write(
                            f"**{codon}** → **{amino_acid}** "
                            f"({full_name}) — **START CODON**"
                        )

                    elif annotation == "Stop codon":

                        st.write(
                            f"**{codon}** → **STOP CODON**"
                        )

                    elif annotation == "Ambiguous codon":

                        st.write(
                            f"**{codon}** → **X** "
                            f"({full_name})"
                        )

                    else:

                        st.write(
                            f"**{codon}** → **{amino_acid}** "
                            f"({full_name})"
                        )

                st.write("### Protein Sequence")

                st.code(protein_sequence)

                st.caption(
                    "X represents an unknown amino acid resulting "
                    "from an ambiguous codon. * represents a stop "
                    "codon in the protein sequence."
                )


        # Tab 4: Compute physical properties like melting temperature and molecular weight
        with tab4:

            st.subheader("Physicochemical Properties")

            st.write("### Melting Temperature")

            tm_col1, tm_col2 = st.columns(2)

            dna_concentration = tm_col1.number_input(
                "DNA concentration (nM)",
                min_value=1.0,
                value=50.0,
                step=1.0
            )

            salt_concentration = tm_col2.number_input(
                "Monovalent salt concentration (mM)",
                min_value=0.1,
                value=50.0,
                step=1.0
            )

            if "N" in sequence:

                st.warning(
                    "Melting temperature cannot be calculated "
                    "because the sequence contains N."
                )

            else:

                tm = calculate_melting_temperature(
                    sequence,
                    dna_concentration,
                    salt_concentration
                )

                st.metric(
                    "Estimated Melting Temperature",
                    f"{tm:.2f} °C"
                )

            st.write("### Molecular Mass")

            if "N" in sequence:

                st.warning(
                    "Molecular mass cannot be calculated because "
                    "the sequence contains N."
                )

            else:

                mass = calculate_molecular_mass(sequence)

                st.metric(
                    "Molecular Mass",
                    f"{mass:.2f} Da"
                )
