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
cd skm-pss-export
pip install .
```
 
For development (editable install):
```bash
git clone https://github.com/NIB-SI/skm-pss-export.git
cd skm-pss-export
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
                     neo4j_password="password",
                     # optional model metadata, written into the exports
                     model_name="PSS", model_version="v3.0.0",
                     creator=["familyName | givenName | organization | email"])

adapter.collect_reactions(access="public")   # connects, queries, closes
adapter.create_sbml(filename="output.sbml")
adapter.create_tabularqual(filename="output.xlsx")
```

## CLI usage

### SBML:

To view the CLI options:
```bash
pss-export to-sbml --help
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

#### SBO terms

The SBO terms are set in `pss_export/pss/pss_export_config.yaml`: per species form (`node_form_to_SBO`),
reaction type (`reaction_subtype_to_SBO`, with the rate law type used for the optional kinetic laws) and
participant role (`node_role_to_SBO`). Choices that differ from what a validator expects:

- **Transporter** (the modifier of a translocation): SBO:0000013 *catalyst*. SBO has no transporter
  participant role; its *transporter* (SBO:0000284) is a physical entity, not a role.
- **Processes** (species forms `process`, `process_active`): SBO:0000375 *process*. SBML expects species
  to be material entities (SBO:0000240), but PSS models processes (e.g. ROS production) as species; the
  term is kept, as it describes them. libSBML warns: `SBO term 'SBO:0000375' on the <species> is not in
  the appropriate branch` (10713).
- **RNAs** (species forms `mrna`, `mirna`, `ncrna`, `ta-sirna`): SBO:0000278 *messenger RNA*,
  SBO:0000316 *microRNA*, SBO:0000334 *non-coding RNA*. These are in SBO's *functional entity* branch
  only, not under *material entity* (unlike *gene*), so libSBML warns as for processes (10713). This is
  an open SBO issue: [EBI-BioModels/SBO#5](https://github.com/EBI-BioModels/SBO/issues/5).

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
