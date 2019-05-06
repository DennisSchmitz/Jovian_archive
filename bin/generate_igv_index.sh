#!/bin/bash

#####################################################################################################################
### To do, ask Robert            ###
###     Example:                                                   ###
#####################################################################################################################

#set -x


#<select>
#  <option value="volvo">Volvo</option>
#  <option value="saab">Saab</option>
#  <option value="mercedes">Mercedes</option>
#  <option value="audi">Audi</option>
#</select>

#9_20689_TCCTGAGCGCGTAAGA_L001_IGVjs.html


htmlfiles="/data/PZN_dev_GitLab/PZN/IGV"

echo "<select>"
for naam in `ls $htmlfiles/*_IGVjs.html`
do
   fn="$(basename ${naam}| awk '{gsub("_IGVjs.html","");print $0}')"
   echo "<option value="http://$(hostname -f)/igv${naam}">${fn}</option>"
done
echo "</select>"
