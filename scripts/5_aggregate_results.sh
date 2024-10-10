!# bin/bash

output_file_name=""

echo "id    chrY_norm" > $output_file_name

cat outfiles/normYcov.*.txt >> output_file_name