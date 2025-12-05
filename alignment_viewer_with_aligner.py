#!/usr/bin/env python3
"""
Browser-based sequence alignment viewer with alignment capability.
"""

import streamlit as st
import tempfile
import os
import subprocess
from pathlib import Path
from collections import defaultdict

def read_aligned_fasta_from_string(fasta_content):
    """Read aligned FASTA from string."""
    sequences = {}
    current_name = None
    current_seq = []
    
    for line in fasta_content.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_name:
                sequences[current_name] = ''.join(current_seq)
            current_name = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_name:
        sequences[current_name] = ''.join(current_seq)
    
    return sequences

def check_if_aligned(sequences):
    """Check if sequences are already aligned (all same length)."""
    if not sequences:
        return False
    lengths = [len(seq) for seq in sequences.values()]
    return len(set(lengths)) == 1

def align_sequences_mafft(fasta_content):
    """Align sequences using MAFFT."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        input_file = f.name
    
    output_file = input_file + '.aligned'
    
    try:
        # Run MAFFT
        result = subprocess.run(
            ['mafft', '--auto', '--thread', '-1', input_file],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode != 0:
            raise Exception(f"MAFFT failed: {result.stderr}")
        
        aligned_content = result.stdout
        
        # Cleanup
        os.unlink(input_file)
        
        return aligned_content, None
        
    except FileNotFoundError:
        return None, "MAFFT not found. Please install MAFFT first."
    except subprocess.TimeoutExpired:
        return None, "Alignment timed out (>5 minutes). Try with fewer/shorter sequences."
    except Exception as e:
        return None, str(e)
    finally:
        if os.path.exists(input_file):
            os.unlink(input_file)

def align_sequences_muscle(fasta_content):
    """Align sequences using MUSCLE."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        input_file = f.name
    
    output_file = input_file + '.aligned'
    
    try:
        # Run MUSCLE (v5 syntax)
        result = subprocess.run(
            ['muscle', '-align', input_file, '-output', output_file],
            capture_output=True,
            text=True,
            timeout=300
        )
        
        if result.returncode != 0:
            # Try MUSCLE v3 syntax
            result = subprocess.run(
                ['muscle', '-in', input_file, '-out', output_file],
                capture_output=True,
                text=True,
                timeout=300
            )
            
            if result.returncode != 0:
                raise Exception(f"MUSCLE failed: {result.stderr}")
        
        with open(output_file, 'r') as f:
            aligned_content = f.read()
        
        # Cleanup
        os.unlink(input_file)
        os.unlink(output_file)
        
        return aligned_content, None
        
    except FileNotFoundError:
        return None, "MUSCLE not found. Please install MUSCLE first."
    except subprocess.TimeoutExpired:
        return None, "Alignment timed out (>5 minutes). Try with fewer/shorter sequences."
    except Exception as e:
        return None, str(e)
    finally:
        if os.path.exists(input_file):
            os.unlink(input_file)
        if os.path.exists(output_file):
            os.unlink(output_file)

def find_differences(sequences):
    """Find positions where sequences differ (ignoring gaps and ambiguity codes)."""
    if not sequences:
        return [], [], []
    
    seq_list = list(sequences.values())
    align_len = len(seq_list[0])
    
    ambiguity_codes = set('RYSWKMBDHVN')
    
    differences = []
    gap_positions = []
    ambiguous_positions = []
    
    for pos in range(align_len):
        bases = [seq[pos].upper() for seq in seq_list]
        bases_set = set(bases)
        
        has_gap = '-' in bases_set
        has_ambiguity = bool(bases_set & ambiguity_codes)
        
        if has_gap:
            gap_positions.append(pos)
        elif has_ambiguity:
            ambiguous_positions.append(pos)
        else:
            if len(bases_set) > 1:
                differences.append(pos)
    
    return differences, gap_positions, ambiguous_positions

def generate_html_visualization(sequences, differences, gap_positions, ambiguous_positions):
    """Generate HTML visualization as string."""
    seq_names = list(sequences.keys())
    align_len = len(list(sequences.values())[0])
    
    html = """
<!DOCTYPE html>
<html>
<head>
    <title>Alignment Viewer</title>
    <style>
        body {
            font-family: 'Courier New', monospace;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .stats {
            background-color: white;
            padding: 15px;
            margin-bottom: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .alignment {
            background-color: white;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .seq-name {
            display: inline-block;
            width: 300px;
            font-weight: bold;
            background-color: #ecf0f1;
            padding: 2px 5px;
            margin: 2px 0;
        }
        .seq-line {
            margin: 2px 0;
            white-space: nowrap;
        }
        .match {
            color: #27ae60;
        }
        .diff {
            background-color: #e74c3c;
            color: white;
            font-weight: bold;
            padding: 0 2px;
        }
        .gap {
            color: #95a5a6;
        }
        .ambig {
            background-color: #f39c12;
            color: white;
            font-weight: bold;
            padding: 0 2px;
        }
        .position-marker {
            color: #7f8c8d;
            font-size: 10px;
            margin-top: 15px;
        }
        .legend {
            color: #7f8c8d;
            margin-bottom: 15px;
        }
    </style>
</head>
<body>
    <div class="stats">
        <h2>📊 Summary Statistics</h2>
"""
    
    html += f"        <p><strong>Number of sequences:</strong> {len(sequences)}</p>\n"
    html += f"        <p><strong>Alignment length:</strong> {align_len:,} bp</p>\n"
    html += f"        <p><strong>Positions with gaps:</strong> {len(gap_positions):,}</p>\n"
    html += f"        <p><strong>Positions with ambiguity codes:</strong> {len(ambiguous_positions):,}</p>\n"
    usable_positions = align_len - len(gap_positions) - len(ambiguous_positions)
    html += f"        <p><strong>Usable positions:</strong> {usable_positions:,}</p>\n"
    html += f"        <p><strong>SNPs (snp-dists method):</strong> {len(differences):,}</p>\n"
    if usable_positions > 0:
        html += f"        <p><strong>Identity (usable positions):</strong> {100 * (1 - len(differences)/usable_positions):.2f}%</p>\n"
    
    html += """
    </div>
    
    <div class="alignment">
        <h2>🧬 Alignment Visualization</h2>
        <div class="legend">
            Legend: <span class="match">Match</span> | 
            <span class="diff">SNP</span> | 
            <span class="gap">Gap</span> | 
            <span class="ambig">Ambiguous</span>
        </div>
"""
    
    diff_set = set(differences)
    gap_set = set(gap_positions)
    ambig_set = set(ambiguous_positions)
    
    chunk_size = 100
    for chunk_start in range(0, align_len, chunk_size):
        chunk_end = min(chunk_start + chunk_size, align_len)
        
        html += f"        <div style='margin: 20px 0;'>\n"
        html += f"        <div class='position-marker'>Position {chunk_start + 1}-{chunk_end}</div>\n"
        
        for name in seq_names:
            seq = sequences[name]
            chunk = seq[chunk_start:chunk_end]
            
            html += f"        <div class='seq-line'>\n"
            html += f"            <span class='seq-name'>{name[:50]}</span> "
            
            for i, base in enumerate(chunk):
                pos = chunk_start + i
                if pos in ambig_set:
                    html += f"<span class='ambig'>{base}</span>"
                elif pos in diff_set:
                    html += f"<span class='diff'>{base}</span>"
                elif base == '-':
                    html += f"<span class='gap'>{base}</span>"
                else:
                    html += f"<span class='match'>{base}</span>"
            
            html += "\n        </div>\n"
        
        html += "        </div>\n"
    
    html += """
    </div>
</body>
</html>
"""
    
    return html

# Streamlit App
st.set_page_config(page_title="Sequence Alignment Viewer", page_icon="🧬", layout="wide")

st.title("🧬 Sequence Alignment Viewer")
st.markdown("Upload sequences to align and visualize differences")

# Sidebar options
with st.sidebar:
    st.header("⚙️ Settings")
    
    alignment_tool = st.selectbox(
        "Alignment Tool",
        ["MAFFT (recommended)", "MUSCLE"],
        help="Choose alignment algorithm. MAFFT is faster and recommended for most cases."
    )
    
    auto_align = st.checkbox(
        "Auto-align if needed",
        value=True,
        help="Automatically align sequences if they're not already aligned"
    )
    
    st.markdown("---")
    st.markdown("""
    **Installation:**
    ```bash
    # MAFFT (recommended)
    conda install -c bioconda mafft
    # or
    sudo apt install mafft
    
    # MUSCLE (alternative)
    conda install -c bioconda muscle
    ```
    """)

# File uploader
uploaded_file = st.file_uploader(
    "Choose a FASTA file",
    type=['fasta', 'fa', 'fna', 'txt'],
    help="Upload aligned or unaligned FASTA sequences"
)

if uploaded_file is not None:
    # Read the file content
    content = uploaded_file.read().decode('utf-8')
    
    # Parse sequences
    try:
        sequences = read_aligned_fasta_from_string(content)
        
        if not sequences:
            st.error("No sequences found in file")
        elif len(sequences) < 2:
            st.error("Need at least 2 sequences to compare")
        else:
            # Check if sequences are aligned
            is_aligned = check_if_aligned(sequences)
            
            if not is_aligned:
                if auto_align:
                    st.info("🔄 Sequences are not aligned. Aligning now...")
                    
                    # Choose alignment tool
                    if "MAFFT" in alignment_tool:
                        aligned_content, error = align_sequences_mafft(content)
                    else:
                        aligned_content, error = align_sequences_muscle(content)
                    
                    if error:
                        st.error(f"Alignment failed: {error}")
                        st.stop()
                    else:
                        st.success("✅ Alignment complete!")
                        sequences = read_aligned_fasta_from_string(aligned_content)
                        
                        # Offer download of aligned sequences
                        st.download_button(
                            label="📥 Download Aligned FASTA",
                            data=aligned_content,
                            file_name="aligned_sequences.fasta",
                            mime="text/plain"
                        )
                else:
                    st.warning("⚠️ Sequences are not aligned. Enable 'Auto-align if needed' in settings to align them.")
                    st.stop()
            else:
                st.success("✅ Sequences are already aligned")
            
            # Find differences
            differences, gap_positions, ambiguous_positions = find_differences(sequences)
            
            # Display summary in sidebar
            with st.sidebar:
                st.header("📊 Quick Stats")
                align_len = len(list(sequences.values())[0])
                usable_positions = align_len - len(gap_positions) - len(ambiguous_positions)
                
                st.metric("Sequences", len(sequences))
                st.metric("Alignment Length", f"{align_len:,} bp")
                st.metric("SNPs", len(differences))
                if usable_positions > 0:
                    identity = 100 * (1 - len(differences)/usable_positions)
                    st.metric("Identity", f"{identity:.2f}%")
            
            # Generate and display HTML
            html_content = generate_html_visualization(
                sequences, differences, gap_positions, ambiguous_positions
            )
            
            # Display in iframe
            st.components.v1.html(html_content, height=800, scrolling=True)
            
            # Download button
            st.download_button(
                label="📥 Download HTML Visualization",
                data=html_content,
                file_name="alignment_visualization.html",
                mime="text/html"
            )
            
    except Exception as e:
        st.error(f"Error processing file: {str(e)}")
        st.exception(e)

else:
    # Show example
    st.info("👆 Upload a FASTA file to get started")
    
    with st.expander("ℹ️ About this tool"):
        st.markdown("""
        This tool can both **align** and **visualize** sequences:
        
        **Features:**
        - 🔄 Auto-align unaligned sequences using MAFFT or MUSCLE
        - 🎨 Visualize differences with color coding
        - 📊 Calculate identity and SNP statistics
        - 📥 Download aligned sequences and HTML visualizations
        
        **Visualization:**
        - **SNPs** are highlighted in red
        - **Gaps** are shown in gray
        - **Ambiguous bases** (IUPAC codes) are shown in orange
        - **Matches** are shown in green
        
        **SNP Counting:**
        - Matches `snp-dists` behavior
        - Excludes positions with gaps
        - Excludes positions with ambiguity codes (RYSWKMBDHVN)
        - Only counts true nucleotide differences
        
        **Supported formats:**
        - Aligned or unaligned FASTA files
        - Multiple sequences for comparison
        """)
    
    with st.expander("📝 Example usage"):
        st.markdown("""
        1. Upload your FASTA file (aligned or unaligned)
        2. If unaligned, sequences will be aligned automatically
        3. View the interactive visualization
        4. Download the results
        
        Works great for:
        - Comparing viral genomes
        - Checking sequence variants
        - Quality control of sequencing data
        - Identifying SNPs and mutations
        """)
