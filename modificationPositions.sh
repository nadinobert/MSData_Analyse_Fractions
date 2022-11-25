#!/usr/bin/env bash
result_id=1502
type=Acetyl
filter="AND (description LIKE '%rdhA%' OR description LIKE '%rdhB%' OR description LIKE '%OmeA%' OR description LIKE '%OmeB%' OR description LIKE '%hupL%' OR description LIKE '%hupS%' OR description LIKE '%hupX%')"
filter=""
sqlite3 ms_data.sqlite "select t.result_id, t.masterModifications FROM peptides t WHERE t.result_id = ${result_id} AND t.modifications LIKE '%${type}%' AND t.confidence = 'High' ${filter}" |\
sed -e 's/^.*\[//g' -e 's/[^0-9;]*//g' | awk -F \; '/.+/{for (i=1;i<=NF;i++){print ($i - 24)}}' | sort -nu