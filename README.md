3 parts to the project you will need to compile:
1. Paper: `Math_4MB3_Acute_Triangle_Paper.Rnw`
2. Presentation: `Math_4MB3_Acute_Triangle_Presentation.Rnw`
3. knitr Supplement: `Math_4MB3_Acute_Triangle_supplement.Rnw`

To begin, please ensure that you do not move around the contents of the folder. You may encounter errors if certain folders are saved elsewhere.

To compile the any of 

1.
Make sure knitr is your default program to weave Rnw files
(see RStudio > Preferences > Sweave, also select the option ‘pdflatex’ in ‘Typeset LaTeX into PDF using:’)

2. 
Open main.Rnw in RStudio

3. 
Either press the button “compile PDF” or run the commands:
	library(knitr)
	knit(“main.Rnw”)

6.
If you pressed the button, the PDF should be open in your default viewer. If you ran the commands, open the resulting .tex file with your preferred LaTeX program and compile the PDF.


To count the total words in the document, run the command `texcount Math_4MB3_Acute_Triangle_Paper.tex` to output a word count.
