# skm-pss-export

Export models from the Plant Stress Signalling (PSS) Neo4j knowledge graph to various file formats. 



## Installation
 
Install directly from GitHub:
```bash
pip install git+https://github.com/NIB-SI/skm-pss-export.git
```
 
Alternatively, clone and install locally:
```bash
git clone https://github.com/NIB-SI/skm-pss-export.git
cd pss-export
pip install .
```
 
For development (editable install):
```bash
git clone https://github.com/NIB-SI/skm-pss-export.git
cd pss-export
pip install -e .
```
 

## PSS database

For testing, a snapshot of the PSS database is available at:
https://github.com/NIB-SI/skm-neo4j

### Connection to Neo4j

You can pass the connection settings (uri, username, password) to the CLI using command-line arguments: 

```bash
  --neo4j-uri TEXT       Neo4j connection URI.
  --neo4j-user TEXT      Neo4j username.
  --neo4j-password TEXT  Neo4j password.
```

Alternatively, set them as environment variables, or in an `.env` file in the current directory, using the following variables:
```bash
MY_NEO4J_URI=bolt://localhost:7687
MY_NEO4J_USER=neo4j
MY_NEO4J_PASSWORD=password
```

Each setting is taken from the argument if given, otherwise from the environment, otherwise from the `.env` file.

If you used the defaults in the [skm-neo4j](https://github.com/NIB-SI/skm-neo4j) repo, you can use the `.env.example` file as provided. 
```bash
mv .env.example .env
```

## Python usage

The package handles the database connection itself: it connects only while collecting data, and closes the connection afterwards.

```python
from pss_export import PSSAdapter

# connection settings as arguments, or from the environment / .env
adapter = PSSAdapter(neo4j_uri="bolt://localhost:7687",
                     neo4j_user="neo4j",
                     neo4j_password="password")

adapter.collect_reactions(access="public")   # connects, queries, closes
adapter.create_sbml(filename="output.sbml")
adapter.create_tabularqual(filename="output.xlsx")
```

## CLI usage

### SBML:

To view the CLI options:
```bash
python pss_adapter_cli.py to-sbml --help
```

Create an SBML file:
```bash
pss-export to-sbml output.sbml --access public
```

Using the model-fixing functions (needs the optional dependencies: `pip install ".[model-fixing]"`):
```bash
pss-export to-sbml output-model-fixes.sbml \
  --access public \
  --model-fixes-identify \
  --model-fixes-apply
```

To add equations to the SBML file, you can use SBMLsqueezer from 
https://github.com/draeger-lab/SBMLsqueezer, e.g.
```bash
java -jar  /path..to..jar/SBMLsqueezer-2.2.jar --sbml-in-file output.sbml  --sbml-out-file output-squeezed
```

#### Known limitations

- **Missing equations**: The SBML file does not contain reaction equations. Add them using SBMLsqueezer as shown above. See also [ThesisRepository](https://github.com/R4d0slav/ThesisRepository).
- **Database consistency**: Continuous updates to the PSS database mean the model may not be fully connected or consistent with modeling assumptions.
- **Disconnected molecules**: Some molecules may be disconnected in the SBML file if:
  - They occur in multiple compartments but lack transport reactions
  - A protein is formed by translation but not activated by an activation reaction
  - A complex is formed but not activated by an activation reaction

### TabularQual:

View available options:
```bash
pss-export to-tabularqual --help
```
 
Create a TabularQual file:
```bash
pss-export to-tabularqual output.xlsx --access public
```

#### Known limitations

- **Database consistency**: Continuous updates to the PSS database mean the model may not be fully connected or consistent with modeling assumptions.
- **Disconnected molecules**: Some molecules may be disconnected in the TabularQual file if:
  - They occur in multiple compartments but lack transport reactions
  - A protein is formed by translation but not activated by an activation reaction
  - A complex is formed but not activated by an activation reaction


## Tests

Two test files are in `tests/`:

- `test_pss_export.py` — unit tests for entity classes, reaction logic, and annotation handling; no database required
- `test_integration.py` — integration tests for SBML and TabularQual export; requires a running Neo4j instance

Run with:
```bash
pytest tests/ -v
```

Integration tests skip automatically if no database connection is available. Set connection details via environment variables or a `.env` file (see [Connection to Neo4j](#connection-to-neo4j)).

> *Tests were AI-generated (Claude) and minimally human-curated.*


## Related Projects
 
- [skm-neo4j](https://github.com/NIB-SI/skm-neo4j) - PSS database snapshots
- [Stress Knowledge Map (SKM)](https://skm.nib.si) - Web interface for the knowledge graph

## License

MIT License. See [LICENSE](LICENSE) for details.
