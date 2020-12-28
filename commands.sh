tshark -r a.pcap -Tfields -e _ws.col.Time -e _ws.col.Length -Y "tcp" > out.txt

tail -n+3 *.csv > mm.csv

awk 'FNR > 1' *.csv > allResults.csv
