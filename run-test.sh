# Test script for UNIX platforms

# Do a complete rebuild because make doesn't always detect when recompilation
# is necessary, and with ccache this is fast anyway.
make clean && make release -j8

# Run
build/sim model_parameters.txt farauti_parameters.txt koliensis_parameters.txt punctulatus_parameters.txt intervention_coverage.txt output.txt 1548338683

echo
echo "Testing for output differences..."
diff -q sample-output/output.txt output.txt && echo "No differences!"
