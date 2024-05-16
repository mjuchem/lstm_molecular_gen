# LSTM Molecular Gen

Trains LSTM with a dataset of SMILES molecular sequences, and uses the model to output new sequences.



Part of the code derived from:
https://github.com/topazape/LSTM_Chem
https://github.com/mattroconnor/deep_learning_coronavirus_cure

# Running

The notebooks can be run under Jupyter Lab.

# Dependencies
- matplotlib
- numpy
- openbabel
- pandas
- rdkit
- tensorflow

# Sandboxed Container

A sandboxed container is provided with the complete environment needed to run the notebooks without installing additional software.

Build the docker image in Linux with `make`, `make build` or `./build.sh`.

Run the docker container with `make run` (ensures image is built beforehand) or `./run.sh`.

Open the lab's URL on your web browser (default: http://localhost:52019).

Configure Junyper Lab with environment variables:
- `NOTEBOOKS_DIR`: notebooks root directory (default: repo root)
- `LAB_ADDR`: host address to bind the server to (default: `127.0.0.1`)
- `LAB_PORT`: host port to bind the server to (default: `52019`)

E.g.:
```
# run the lab on the default address and port
make run

# run the lab on http://127.0.0.2:51776
make run LAB_ADDR=127.0.0.2 LAB_PORT=51776
```
