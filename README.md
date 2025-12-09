# Sequence Alignment Viewer

A browser-based tool for aligning and visualizing DNA/RNA sequence alignments with automatic SNP detection and beautiful color-coded output.

## Features

- **Auto-align** unaligned FASTA sequences using MAFFT or MUSCLE
- **Color-coded visualization** of matches, SNPs, gaps, and ambiguous bases
- **Detailed statistics** including identity percentage and SNP counts
- **Download results** as HTML visualizations or aligned FASTA files
- **100% local** - all processing happens on your machine, no data uploaded
- **Easy to use** - simple drag-and-drop interface

## What It Does

This tool helps you:
1. Compare multiple DNA/RNA sequences
2. Identify single nucleotide polymorphisms (SNPs)
3. Visualize alignment gaps and ambiguous bases
4. Calculate sequence identity metrics
5. Export publication-ready visualizations

Perfect for viral genomics, outbreak investigation, and comparative sequence analysis.

## Prerequisites

### Required Software

1. **Python 3.8 or higher**
   - Check your version: `python3 --version`
   - Download from: https://www.python.org/downloads/

2. **pip** (Python package installer)
   - Usually comes with Python
   - Check: `pip --version` or `pip3 --version`

### Installation Steps

#### Step 1: Install Python Dependencies

Open a terminal/command prompt and run:

```bash
pip install streamlit
```

Or if you're using `pip3`:

```bash
pip3 install streamlit
```

#### Step 2: Install Alignment Software (Optional but Recommended)

If you want to align unaligned sequences, install MAFFT:

**Option A: Using Conda (Recommended)**
```bash
# If you have conda/mamba installed
conda install -c bioconda mafft

# Or using mamba (faster)
mamba install -c bioconda mafft
```

**Option B: Using Package Manager**

On **Ubuntu/Debian**:
```bash
sudo apt update
sudo apt install mafft
```

On **macOS** (using Homebrew):
```bash
brew install mafft
```

On **Windows**:
- Download from: https://mafft.cbrc.jp/alignment/software/windows.html
- Or use WSL (Windows Subsystem for Linux) and follow Ubuntu instructions

**Alternative: MUSCLE**
```bash
# Using conda
conda install -c bioconda muscle

# Or download from: http://www.drive5.com/muscle/
```

**To verify installation:**
```bash
mafft --version
# Should show version number if installed correctly
```

**Note:** If you only want to visualize already-aligned sequences, you can skip the alignment software installation.

#### Step 3: Download the Tool

Download `alignment_viewer_script.py` and save it to a folder on your computer.

## Usage

### Starting the Tool

1. Open a terminal/command prompt
2. Navigate to the folder containing the script:
   ```bash
   cd /path/to/folder
   ```
3. Run the application:
   ```bash
   streamlit run alignment_viewer_script.py
   ```
4. Your browser will automatically open to `http://localhost:8501`

### Using the Tool

#### For Aligned Sequences:
1. Click "Browse files" or drag your aligned FASTA file
2. View the visualization immediately
3. Download HTML or aligned FASTA if needed

#### For Unaligned Sequences:
1. Make sure "Auto-align if needed" is checked (default: ON)
2. Upload your unaligned FASTA file
3. The tool will automatically align your sequences
4. View the results and download if needed

#### Settings (in sidebar):
- **Alignment Tool**: Choose between MAFFT (faster) or MUSCLE
- **Auto-align if needed**: Toggle automatic alignment on/off

### Example Input File

Your FASTA file should look like this:

```
>Sequence_1
ATGGCTAGCTAGCTAGCTAG
>Sequence_2
ATGGCTAGCTGGCTAGCTAG
>Sequence_3
ATGGCTAGCTAGCTAGCGAG
```

Sequences can be:
- Unaligned (different lengths) - will be aligned automatically
- Already aligned (same lengths) - will be visualized directly
- DNA or RNA
- Any number of sequences (2 or more)

## Understanding the Output

### Color Coding:
- **Green (Match)**: Position matches across all sequences
- **Red (SNP)**: Single nucleotide polymorphism - sequences differ here
- **Gray (Gap)**: Insertion/deletion gap (-)
- **Orange (Ambiguous)**: IUPAC ambiguity code (N, R, Y, etc.)

### Statistics Explained:
- **Alignment length**: Total length of the aligned sequences
- **Positions with gaps**: Number of positions containing gaps
- **Positions with ambiguity codes**: Positions with N, R, Y, etc.
- **Usable positions**: Positions without gaps or ambiguity (used for SNP counting)
- **SNPs**: Number of true nucleotide differences (snp-dists compatible)
- **Identity**: Percentage of usable positions that match

### SNP Counting Method:
This tool uses the same method as `snp-dists`:
- Only counts positions with standard bases (A, T, G, C)
- Excludes positions with gaps (-)
- Excludes positions with ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N)

## Troubleshooting

### "ModuleNotFoundError: No module named 'streamlit'"
**Solution:** Install streamlit:
```bash
pip install streamlit
```

### "MAFFT not found" or "MUSCLE not found"
**Solution:** Either:
- Install the alignment tool (see Step 2 above), or
- Use only pre-aligned sequences (no installation needed)

### "Port 8501 is already in use"
**Solution:** Either:
- Close the other streamlit app, or
- Use a different port:
  ```bash
  streamlit run alignment_viewer_script.py --server.port 8502
  ```

### Browser doesn't open automatically
**Solution:** Manually open your browser and go to:
```
http://localhost:8501
```

### Permission denied errors (Linux/Mac)
**Solution:** You may need to use `sudo` for system-wide installations:
```bash
sudo apt install mafft
```
Or use conda which doesn't require sudo.

### Alignment takes too long
**Possible causes:**
- Very large sequences (>100kb each)
- Many sequences (>50)

**Solutions:**
- Use MAFFT instead of MUSCLE (faster)
- Pre-align sequences using command-line tools
- Break into smaller batches

## Privacy & Security

- **100% local processing** - no data leaves your computer
- **No internet required** after installation
- **Files processed in memory** - not saved to disk
- **Works offline** once dependencies are installed
- **Safe for sensitive/private data**

Perfect for working with:
- Patient sequences
- Unpublished research data
- Proprietary genomic data
- Data with privacy restrictions

## Advanced Usage

### Custom Port
```bash
streamlit run alignment_viewer_script.py --server.port 8080
```

### Headless Mode (server without browser)
```bash
streamlit run alignment_viewer_script.py --server.headless true
```

### Command Line Alignment (without GUI)
If you want to use the original command-line script:
```bash
python view_alignment.py input.fasta output.html
```

## System Requirements

- **RAM**: 2GB minimum, 4GB+ recommended for large alignments
- **Processor**: Any modern CPU
- **Disk Space**: ~500MB for Python + dependencies
- **Operating System**: Windows, macOS, or Linux

## File Size Limits

Recommended maximums for good performance:
- **Sequence length**: <100,000 bp each
- **Number of sequences**: <100 sequences
- **Total file size**: <50MB

Larger datasets will work but may take longer to process.

## Support

### Common Issues

**Q: Can I use this for protein sequences?**
A: Not currently - it's designed for DNA/RNA. Modifications would be needed for amino acids.

**Q: Can multiple people use this at once?**
A: Each person needs to run their own instance. It's single-user by design.

**Q: Does this work with BAM/VCF files?**
A: No, only FASTA format is supported.

**Q: Can I deploy this on a server for my team?**
A: Yes! You can run it on an internal server and share the URL with your team (within your network).

### Getting Help

If you encounter issues:
1. Check the Troubleshooting section above
2. Verify all prerequisites are installed
3. Check that your FASTA file is properly formatted
4. Try with a small test file first

## Citation

If you use this tool in your research, please cite the underlying tools:

**For MAFFT alignments:**
```
Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment 
software version 7: improvements in performance and usability. 
Molecular Biology and Evolution, 30(4), 772-780.
```

**For MUSCLE alignments:**
```
Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high 
accuracy and high throughput. Nucleic Acids Research, 32(5), 1792-1797.
```

## License

This tool is provided for research and educational use.

## Version History

- **v1.1** - Added auto-alignment capability with MAFFT/MUSCLE
- **v1.0** - Initial release with visualization features

## Tips for Best Results

1. **Sequence naming**: Use descriptive names in FASTA headers (avoid special characters)
2. **File preparation**: Remove any non-sequence characters before upload
3. **Quality check**: Remove low-quality regions from sequences before alignment
4. **Batch size**: For >50 sequences, consider splitting into smaller groups
5. **Reference sequence**: Put your reference sequence first in the FASTA file

## Acknowledgments

This tool uses:
- Streamlit for the web interface
- MAFFT/MUSCLE for sequence alignment
- Standard Python libraries for processing

Designed for public health genomic surveillance and molecular epidemiology.

---

**Questions or suggestions?** Feel free to provide feedback!
