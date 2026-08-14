import streamlit as st

st.title("seqMetrics: DNA Sequence Analyzer")
st.write("Calculate GC content percentage and length from a DNA sequence string.")

# Input box for DNA sequence..
dna_seq = st.text_area("Enter DNA sequence (5' to 3'):", "").upper().strip()

if dna_seq:
    # Clean sequence in case of whitespace or numbers....
    valid_bases = set("ATGC")
    if all(base in valid_bases for base in dna_seq):
        seq_length = len(dna_seq)
        g_count = dna_seq.count('G')
        c_count = dna_seq.count('C')
        gc_content = ((g_count + c_count) / seq_length) * 100 if seq_length > 0 else 0
        
        st.subheader("Results")
        st.metric(label="Sequence Length", value=f"{seq_length} bp")
        st.metric(label="GC Content", value=f"{gc_content:.2f}%")
    else:
        st.error("Invalid DNA sequence. Please use only standard nucleotide characters: A, T, G, C.")
      
