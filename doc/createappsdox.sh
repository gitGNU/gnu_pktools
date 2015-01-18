#!/bin/bash

echo "create header"

echo "create dox files for applications"
for file in ../src/apps/pk*.cc;do    
    THETOOL=$(basename $file .cc)
    echo ${THETOOL}

    THESHORTDESCRIPTION=$(grep "${THETOOL}.cc: " $file | awk -v FS=':' '{print $2}')
    USAGE=$(${THETOOL} -h|grep Usage)
    cat > ../doc/${THETOOL}.dox <<EOF
\section $THETOOL $THETOOL
$THESHORTDESCRIPTION
## SYNOPSIS

<code>
  ${USAGE}
</code>

EOF

    #create description if not exists
    if [ ! -f ../doc/description_${THETOOL}.dox ];then
	touch ../doc/description_${THETOOL}.dox
    fi
    cat ../doc/description_${THETOOL}.dox >> ../doc/${THETOOL}.dox

## OPTIONS ##
    cat >> ../doc/${THETOOL}.dox <<EOF
\section ${THETOOL}_options Options
 - use either \`-short\` or \`--long\` options (both \`--long=value\` and \`--long value\` are supported)
 - short option \`-h\` shows basic options only, long option \`--help\` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
EOF

    ${THETOOL} --doxygen|sed '$d' >> ../doc/${THETOOL}.dox
    echo >> ../doc/${THETOOL}.dox
    if [ -f examples_${THETOOL}.dox ];then
	echo "Examples" >> ../doc/${THETOOL}.dox
	echo "========" >> ../doc/${THETOOL}.dox
	echo "Some examples how to use ${THETOOL} can be found \\ref examples_${THETOOL} \"here\"" >> ../doc/${THETOOL}.dox
    fi
    if [ -f faq_${THETOOL}.dox ];then
	echo "FAQ" >> ../doc/${THETOOL}.dox
	echo "========" >> ../doc/${THETOOL}.dox
	echo "Frequently asked questions on ${THETOOL} can be found \\ref faq_${THETOOL} \"here\"" >> ../doc/${THETOOL}.dox
    fi
done

echo "create general dox file for aps list"
echo "\section available_tools Available tools" > ../doc/apps.dox
for file in ../src/apps/pk*.cc;do
    THETOOL=$(basename $file .cc)
    THESHORTDESCRIPTION=$(grep "${THETOOL}.cc: " $file | awk -v FS=':' '{print $2}')
    echo "- \\ref ${THETOOL} ${THESHORTDESCRIPTION}"; 
done >> ../doc/apps.dox

#remove depricated utilities and those not ready to publish"

for TOOL in pkeditogr pkenhance pkkalman pkndvi pkreclass; do 
    rm -f ../doc/${TOOL}.dox ../html/md_doc_${TOOL}.html ../html/${TOOL}_8cc_source.html
    sed -i "/${TOOL}/d" ../doc/apps.dox
done
echo "Savannah repository for homepage can only be maintained via cvs"
#mkdir ~/tmp
#cd ~/tmp
#cvs -z3 -d:ext:kempenep@cvs.sv.gnu.org:/web/pktools co pktools
#cd pktools/html
#rm *
#cvs rm *
#rsync -avz <orig html>/ ~/tmp/pktools/html
#cvs add *.*
#cvs commit -m "update of repository homepage"
#rm -r ~/tmp/pktools
echo "ftp to downloads"
#sftp kempenep@download.savannah.gnu.org:/releases/pktools
