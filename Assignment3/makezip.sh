result=${PWD##*/}
jupyter nbconvert $result.ipynb --to latex --template clean_report.tplx
pdflatex $result.tex
rm $result.zip
zip $result.zip $result.ipynb $result.tex $result.pdf ${PWD##*/}_files
