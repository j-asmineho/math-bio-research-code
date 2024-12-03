3 parts to the project you will need to compile:
1. Paper
2. Presentation
3. knitr Supplement

To begin, please ensure that you do not move around the contents of the folder. You may encounter errors if certain folders are saved elsewhere.

To compile the paper
1. Ensure that the Images folder has been downloaded and is located in the same place as Math_4MB3_Acute_Triangle_Paper.tex
2. Press "Typset" in LaTex or run the following commands:

     render pdflatex Math_4MB3_Acute_Triangle_Paper.tex

To compile the knitr Supplement
1. Make sure knitr is your default program to weave Rnw files
(see RStudio > Preferences > Sweave, also select the option ‘pdflatex’ in ‘Typeset LaTeX into PDF using:’)

2. Open main.Rnw in RStudio

3. Either press the button “compile PDF” or run the commands:
	  library(knitr)
	  knit(“supplement.Rnw”)
   
4.If you pressed the button, the PDF should be open in your default viewer. If you ran the commands, open the resulting .tex file with your preferred LaTeX program and compile the PDF.

Presentation
1. Either press "Typeset" button in LaTex or run the command:

     render pdflatex Math_4MB3_Acute_Triangle_Paper.tex

To count the total words in the document, run the command `texcount Math_4MB3_Acute_Triangle_Paper.tex` to output a word count.
