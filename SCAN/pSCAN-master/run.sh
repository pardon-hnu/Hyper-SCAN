declare -a epision=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
# declare -a miu=(2 3 4 5 10 15)
declare -a miu=(6 7 8 9)

for mi in "${miu[@]}"
do
    for ep in "${epision[@]}"
    do
        ./pscan /home/hnu/Disk0/ParDon/hypergraph/Hyper2General/coauth-DBLP $ep $mi output
        echo "------------"
    done
done
