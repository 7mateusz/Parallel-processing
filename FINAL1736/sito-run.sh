#!/usr/bin/env bash
[ $# -ne 4 ] && { echo "usage: $0 <n> <k> <mod> <workers>" >&2; exit 1; }
n=$1; k=$2; mod=$3; w=$4
mkdir -p slices
base="slices/${n}_${k}_${mod}"
touch "$base.done"

for r in $(seq 0 $((mod-1))); do
  grep -qx "$r" "$base.done" && continue
  (
    geng -c $n $k:$k $r/$mod 2>>"${base}_r${r}.log" \
      | ./build/sito-seq > "${base}_r${r}.g6" 2>>"${base}_r${r}.log"
    grep -q '^>Z' "${base}_r${r}.log" && grep -q 'done t=' "${base}_r${r}.log" \
      && echo "$r" >> "$base.done"
  ) &
  while [ "$(jobs -rp | wc -l)" -ge "$w" ]; do wait -n; done
done
wait
echo "$(wc -l < $base.done)/$mod slices done: $(cat ${base}_r*.g6 2>/dev/null | wc -l) survivors, logs in slices/"
