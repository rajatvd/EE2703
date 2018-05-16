result=${PWD##*/}
rm -r ./${PWD##*/}_files
jupyter nbconvert $result.ipynb --to latex --template clean_report.tplx
pdflatex $result.tex
rm $result.zip
zip -r $result.zip $result.ipynb $result.tex $result.pdf *.py ./${PWD##*/}_files
