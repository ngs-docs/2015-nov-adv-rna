trackDb="$1"
lib_assemblies="$2"

lib_assembly_array=()
## add the new data
while read assembly; do
  echo $assembly
  ## ensure the assembly was not loaded before
  if [ $(echo ${lib_assembly_array[@]} | grep -o $assembly | wc -l) -eq 0 ]; then
    lib_assembly_array+=($assembly)
fi; done < $lib_assemblies

## restart the trackDb
> $trackDb

## edit the trackDb
priority=1
for assembly in ${lib_assembly_array[@]}; do
  ## create non-composite entry of the assembly
  echo $assembly;
  filename=$(basename $assembly)
  echo "track $assembly" >> $trackDb
  echo "bigDataUrl BigBed/$filename.BigBed" >> $trackDb
  echo "shortLabel $assembly" >> $trackDb
  echo "longLabel $assembly" >> $trackDb
  echo "type bigBed 12" >> $trackDb
  echo "colorByStrand 255,0,0 0,0,255" >> $trackDb
  echo "visibility dense" >> $trackDb
  echo "priority $priority" >> $trackDb
  echo "html $filename" >> $trackDb
  echo " " >> $trackDb
  let "priority++"
done



