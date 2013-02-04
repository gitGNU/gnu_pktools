#!/bin/bash

echo "create header"
cat > ../doc/header.dox <<EOF
\section thetool thetool
theshortdescription

## SYNOPSIS
<code>
thesynopsis
</code>

## DESCRIPTION ##
thelongdescription

## OPTIONS ##
 - use either \`-short\` or \`--long\` options (both \`--long=value\` and \`--long value\` are supported)
 - short option \`-h\` shows basic options only, long option \`--h\` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
EOF

echo "create dox files for applications"
for file in ../src/apps/pk*.cc;do    
    THETOOL=$(basename $file .cc)
    echo ${THETOOL}
    cat ../doc/header.dox > ../doc/${THETOOL}.dox
    sed -i "s/thetool/$THETOOL/g" ../doc/${THETOOL}.dox
    THESHORTDESCRIPTION=$(grep "${THETOOL}.cc: " $file | awk -v FS=':' '{print $2}')
    sed -i "s/theshortdescription/$THESHORTDESCRIPTION/" ../doc/${THETOOL}.dox
    THESYNOPSIS="${THETOOL} [OPTIONS]"
    sed -i "s/thesynopsis/$THESYNOPSIS/" ../doc/${THETOOL}.dox
    THELONGDESCRIPTION="$THESHORTDESCRIPTION more..."
    sed -i "s/thelongdescription/$THELONGDESCRIPTION/" ../doc/${THETOOL}.dox
    ${THETOOL} --doxygen|sed '$d' >> ../doc/${THETOOL}.dox
done

echo "create general dox file for aps list"
for file in ../src/apps/pk*.cc;do
    THETOOL=$(basename $file .cc)
    THESHORTDESCRIPTION=$(grep "${THETOOL}.cc: " $file | awk -v FS=':' '{print $2}')
    echo "- \\ref ${THETOOL} ${THESHORTDESCRIPTION}"; 
done > ../doc/apps.dox

echo "Savannah repository for homepage can only be maintained via cvs"
#mkdir ~/tmp
#cd ~/tmp
#cvs -z3 -d:ext:kempenep@cvs.sv.gnu.org:/web/pktools co pktools"
#cd pktools/html
#rm *
#cvs rm *
#rsync -avz <orig html>/ ~/tmp/pktools/html
#cvs add *
#cvs commit -m "update of repository homepage"
#rm -r ~/tmp/pktools