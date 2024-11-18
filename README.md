# <p align="center"> Nice to meet you！ :kissing_closed_eyes: </p>
--------
#### This is ****AbnumPro**** , Before You Begin:
***Important note***:Since antibody numbering requires relatively low computational resources and benefits from real-time visualization, most users typically access online antibody numbering platforms on Windows systems. Additionally, the [Cao Lab](http://cao.labshare.cn/AbRSA/) has developed an antibody numbering software for Linux, enabling Linux and Mac users to perform both small- and large-scale numbering tasks, even though it may also function as an online service. Therefore, AbnumPro is currently developed exclusively for ***Windows*** users. Users can download the latest version of AbnumPro for offline use.<br>

***Current version: v1.1.0*** <br> 

If ***Chinese*** is your native language, congratulations~~~🎉🎉🎉 You can use AbnumPro seamlessly! If not, no worries—you can switch the language to English.👻

The AbnumPro is free for anyone to use and distribute, including for commercial and academic purposes. For details, please refer to the License.

This document contains ***no code or scripts***, so feel free to use it with confidence!

----------

## Overview
## Getting Started
### Prerequisites

### Installation
* STEP1: Download the compressed file named `AbnumPro.zip` from the repository.<br>
* STEP2: After extracting the files, no installation is required—simply double-click the `.exe` file to run the program.<br>

### Usage
*  The default language of AbnumPro is Chinese, but you can switch to English from the top-right corner of the software.<br>
  * The left sidebar serves as the AbnumPro function panel, integrating three features: antibody numbering, CDR prediction, and ABR prediction.<br>
    * antibody numbering
        * After inputting the antibody sequence, select the antibody numbering scheme and click the Submit button to perform antibody numbering. You can also click the Example button to populate the input box with sample data.
    * CDR prediction
        * The method is the same as above. Due to scheme limitations, CDR prediction is only available for the Kabat, Chothia, and IMGT schemes.
    * ABR prediction
        * ABR prediction in AbnumPro is based on a Hidden Markov Model and does not require selecting a numbering scheme.
  * All prediction results can be saved by clicking the Save button.

## Performance Evaluation
<div align="center">

|Scheme|Light Chain|Heavy chain|Complete Antibody|
| ---------- | -----------| -----------| -----------|
|  `Kabat`|0.852/0.356|0844/0.482|0.846/0.424|
| `Chothia`   |0.795/0.334| 0.816/0.632   | 0.808/0.473   |
| `IMGT`   |0.744/0.589|0.706/0.594| 0.720/0.592   |
| `ABR`   |0.947/0.343|0.903/0.523|0.920/0.433|

</div>

## License
* This project is licensed under the MIT License - see the [LICENSE]() file for details.
----------
This project is jointly led by [HLAB](https://i.uestc.edu.cn/hlab/Journal.html) and Xie Lab. For any questions, don't hesitate to get in touch with:<br>
* Wenzhen Li, first author, at bioinfowenz@gmail.com.<br>
* Prof. Jian Huang, corresponding author, at hj@uestc.edu.cn.<br>

Please cite:

# <p align="center"> Have a goooood day！ :kissing_closed_eyes: </p>
