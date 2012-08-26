SCRIPTNAME="default_post_thermalization.script.lisp"

for i in {1..28}; do 
    nice sbcl --dynamic-space-size 2000 --script $SCRIPTNAME >> "120824_ensemble_building_T016_$i.log" &
    sleep 2s
    echo $i
done

echo "All finished"
