# TH-GEM Foil Analyzer

**Optical Assessment Tool for Gas Electron Multiplier Detectors**

---

## Overview

FoilAnalyzer is a Python-based software package for automated optical characterization of Thick Gas Electron Multiplier (TH-GEM) foils through high-resolution image analysis. The tool implements color-based segmentation using a nearest-neighbor algorithm to identify and measure internal hole diameters across the detector surface.

This software was developed as part of a quality assurance protocol for TH-GEM foils manufactured through CERN's Knowledge Transfer program, supporting the qualification of detectors for high energy physics experiments.

## Technical Specifications

### Performance Metrics
- **Identification Rate**: over 95%
- **Measurement Precision**: (0.31 ± 0.03) mm for internal diameter
- **Input Resolution**: 3600 DPI recommended
- **Image Format**: 24-bit color JPG

### Algorithm Overview
The analysis pipeline consists of three main stages:

1. **Color-Space Segmentation**: Pixels are classified into copper substrate, external rim, or internal hole regions using Euclidean distance minimization in RGB space
2. **Geometric Analysis**: Vertical column scanning with pitch-based positioning to locate hole centers
3. **Diameter Measurement**: Sequential pixel counting with DPI-based conversion to physical dimensions

## Requirements

### Core Dependencies
```
Python >= 3.7
numpy >= 1.19.0
matplotlib >= 3.3.0
pandas >= 1.1.0
```

### Optional Dependencies
For histogram visualization (Step 3):
```
ROOT >= 6.0 (PyROOT)
```

Installation of PyROOT varies by platform. Consult the [ROOT Installation Guide](https://root.cern/install/) for system-specific instructions.

## Installation

Clone the repository:
```bash
git clone https://github.com/g-lickel/FoilAnalyzer.git
cd FoilAnalyzer
```

Install Python dependencies:
```bash
pip install numpy matplotlib pandas
```

## Usage

### Configuration

Edit the initial parameters in `FoilAnalyzer.py` to match your detector specifications:

```python
# Detector specifications
file = 'CMGL006_10X10_3600_Cor_24Bits'   # Input filename (without extension)
inDiam = 0.3      # Expected internal diameter (mm)
outDiam = 0.5     # Expected external diameter (mm)
pitch = 1.55      # Hole pitch (mm)
dpi = 3600        # Image resolution (dots per inch)
firstCenter = 2.75  # First column center position (mm from edge)
```

### Reference Color Calibration

The `data()` function contains RGB reference values for material classification. These should be calibrated for each scanning setup:

```python
def data():
    copper = [[R, G, B], ...]    # Copper substrate samples
    outhole = [[R, G, B], ...]   # External hole rim samples
    inhole = [[R, G, B], ...]    # Internal hole samples
    return (copper, outhole, inhole)
```

### Execution

Run the main script:
```bash
python3 FoilAnalyzer.py
```

The program provides three operational modes:

#### Option 1: Image Analysis
Processes the raw scan and generates a grayscale-classified image:
- **Input**: `[filename].jpg` (raw scan at specified DPI)
- **Output**: `[filename]_GS.png` (classified image: black=copper, gray=outer rim, white=inner hole)
- **Processing Time**: Varies with image dimensions (~hours for 3600 DPI full foil scans)

#### Option 2: Diameter Measurement
Performs geometric analysis on the classified image:
- **Input**: `[filename]_GS.png` (from Option 1)
- **Output**: `[filename]_hist.csv` (measured diameters in mm, one per line)
- **Algorithm**: Column-wise scanning with pitch-based grid positioning

#### Option 3: Statistical Analysis
Generates histogram with Gaussian fit (requires PyROOT):
- **Input**: `[filename]_hist.csv` (from Option 2)
- **Output**: Interactive ROOT canvas with fitted distribution
- **Parameters**: Mean and standard deviation of hole diameter distribution

### Typical Workflow

```bash
# Step 1: Process raw image
python FoilAnalyzer.py
> 1  # Select image analysis

# Step 2: Measure diameters
python FoilAnalyzer.py
> 2  # Select measurement

# Step 3: Visualize results
python FoilAnalyzer.py
> 3  # Select histogram plotting
```

## Methodology

### Image Acquisition Protocol
- **Scanner**: Flatbed scanner with optical resolution ≥ 3600 DPI
- **Illumination**: Uniform backlight or transmitted light configuration
- **Format**: 24-bit RGB color, JPEG compression level ≥ 90%
- **Calibration**: Include ruler or reference scale in scan

### Color Classification Algorithm

For each pixel at position (x, y):

1. Compute Euclidean distance in RGB space to each reference color set:
   ```
   d_i = sqrt((R - R_ref,i)² + (G - G_ref,i)² + (B - B_ref,i)²)
   ```

2. Assign pixel to class with minimum distance:
   ```
   class(x,y) = argmin_i(d_i)
   ```

3. Map to output values:
   - Copper: [0, 0, 0] (black)
   - External rim: [189, 189, 189] (gray)
   - Internal hole: [255, 255, 255] (white)

### Diameter Extraction Algorithm

Column-based scanning approach:

1. **Grid Definition**: Calculate column centers based on `pitch` and `firstCenter`
2. **Vertical Scanning**: For each column, scan vertically pixel-by-pixel
3. **Hole Detection**: Identify continuous white pixel sequences
4. **Diameter Calculation**: 
   - Count maximum consecutive white pixels per hole
   - Convert to physical units: `diameter_mm = pixel_count × 25.4 / DPI`
5. **Validation**: Filter measurements based on expected `inDiam` range

## Output Files

### Intermediate Results
- `[filename]_GS.png`: Classified image for visual inspection and measurement input

### Analysis Results
- `[filename]_hist.csv`: Raw diameter measurements
  ```
  0.298
  0.312
  0.305
  ...
  ```

### Statistical Summary
Generated by ROOT histogram (Option 3):
- Mean diameter and uncertainty
- Standard deviation
- Distribution histogram with Gaussian fit overlay

## Known Limitations

### Image Acquisition
- **Reflection Artifacts**: Specular reflections from copper surface can cause misclassification
- **Non-uniform Illumination**: Intensity gradients affect color-based segmentation
- **Focus Depth**: Thick substrates may exhibit focus variation across scan area

### Algorithm Constraints
- **Fixed Grid Assumption**: Requires regular hole pattern with known pitch
- **Color Dependency**: Sensitive to scanner calibration and color profile
- **Edge Effects**: Holes near foil boundaries may be partially captured

### Computational Performance
- Memory usage scales with image dimensions (~16 GB for 3600 DPI, 10×10 cm² scan)
- Processing time scales as O(n²) with pixel count

## Future Development

### Planned Enhancements
- Machine learning-based classification (KNN architecture)
- Adaptive threshold determination
- Pitch variation accommodation
- Parallel processing implementation
- Graphical user interface
- Real-time processing pipeline

### Research Directions
- Multi-scale analysis for hierarchical detector structures
- 3D reconstruction from focus stacking
- Automated defect classification
- Integration with electrical characterization data

## Validation and Quality Metrics

### Accuracy Assessment
The tool's performance has been validated against:
- Manual measurements using calibrated microscopy
- Manufacturer specifications for THGEM geometry
- Statistical consistency across multiple scans of identical foils

### Quality Control Parameters
- Identification rate > 95% considered acceptable
- Diameter distribution width < 0.05 mm indicates good measurement precision
- Visual inspection of classified image required to verify segmentation quality

## Citation

If you use FoilAnalyzer in your research, please cite:

```bibtex
@software{lickel2025foilanalyzer,
  author = {Lickel, Guilherme V. and Manganote, Edmilson J. T.},
  title = {FoilAnalyzer: Optical Assessment Tool for TH-GEM Foils},
  year = {2025},
  url = {https://github.com/g-lickel/FoilAnalyzer},
  institution = {Instituto de Física Gleb Wataghin, Universidade Estadual de Campinas}
}
```

## References

1. Colaleo, A., et al. (2015). *CMS Technical Design Report for the Muon Endcap GEM Upgrade*. CERN-LHCC-2015-012, CMS-TDR-013.

2. Sauli, F. (2016). *The gas electron multiplier (GEM): Operating principles and applications*. Nuclear Instruments and Methods in Physics Research Section A, 805:2–24. doi:[10.1016/j.nima.2015.07.060](https://doi.org/10.1016/j.nima.2015.07.060)

## Acknowledgments

This project was developed at the Cosmic Rays and Chronology Department (DRCC), Gleb Wataghin Institute of Physics, UNICAMP.

**Funding**: This study was financed, in part, by the São Paulo Research Foundation (FAPESP), Brasil. Process Number #2023/11934-6.

**Contributors**:
- Guilherme V. Lickel (Development, Algorithm Design)
- Prof. Edmilson J. T. Manganote, PhD (Supervision, Methodology)

## Contact

**Guilherme V. Lickel**  
Gleb Wataghin Institute of Physics
Universidade Estadual de Campinas  
Email: glickel@ifi.unicamp.br

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Repository**: [github.com/g-lickel/FoilAnalyzer](https://github.com/g-lickel/FoilAnalyzer)  
**Issues**: Report bugs or request features via GitHub Issues  
**Version**: 1.0.0 (March 2025)
