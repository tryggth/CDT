for i in {1..24}; do 
    nohup nice sbcl --dynamic-space-size 2000 --script $1 >> $1.log &
    sleep 2s
    echo $i
done

echo "All finished"
