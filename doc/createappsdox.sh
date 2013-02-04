#!/bin/bash

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

for file in ../src/apps/pk*.cc;do
    THETOOL=$(basename $file .cc)
    THESHORTDESCRIPTION=$(grep "${THETOOL}.cc: " $file | awk -v FS=':' '{print $2}')
    echo "- \\ref ${THETOOL} ${THESHORTDESCRIPTION}"; 
done > ../doc/apps.dox