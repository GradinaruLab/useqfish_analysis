<div id="top"></div>
<!-- PROJECT LOGO -->
<br />
<div align="center">
<!--   <a href="https://github.com/minjeejang/useqfish_analysis">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a> -->

  <h3 align="center">USeqFISH Analysis</h3>

  <p align="center">
    Ultrasensitive, sequential flourescence in situ hybridization (USeqFISH) is a new, highly sensitive spacial transcriptomics method.
    This is data analysis pipeline for images generated in USeqFISH.  
    <br />
    <a href="https://github.com/minjeejang/useqfish_analysis"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    ·
    <a href="https://github.com/minjeejang/useqfish_analysis/issues">Report Bug</a>
    ·
    <a href="https://github.com/minjeejang/useqfish_analysis/issues">Request Feature</a>
  </p>
</div>


### Built With

* [Python](https://www.python.org/)
* [Cellpose](https://www.cellpose.org/)
* [Dask](https://dask.org/)

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.


### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/minjeejang/useqfish_analysis.git
   ```
3. Install Requirements
   ```sh
   pip install --upgrade pip
   pip install -r requirements.txt 
   ```
   
<!-- USAGE EXAMPLES -->
## Usage

To run the pipeline, you'll need multiple round images in .tif format.

#### 1. Set the parameters.
Open ```params.py``` and set the parameters for your project.
* nC: Number of channels each round has
* nR: Number of total rounds
* roundRef: Which round contains the cell reference
* cellCh: channel containing cell reference
* image_size: image size in pixels
* sigma: estimate size of spots to detect
* color_shifts: channel shifts for each channel, provided by the equipment

#### 2. Run the pipeline
 ```sh
   python get_expression_matrix.py [path]
   ```
   [path] is the path to the directory containting your images. All .tif images in this path will be analyzed.
   Note: script will take a while (30 minutes+) to run, depending on image size and machine specifications.
  
#### 3.  Examine results
  Once the run is finished, there will be a ```results.xlsx``` file in your image path. This contains a the expression by cell ID.
  
<p align="right">(<a href="#top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

Min Jee Jang

Project Link: [https://github.com/minjeejang/useqfish_analysis/](https://github.com/minjeejang/useqfish_analysis/)

<p align="right">(<a href="#top">back to top</a>)</p>
