import yaml
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Tuple, List, Optional

IDENTIFIERS_ORG = "http://identifiers.org"


@dataclass(frozen=True)
class Annotation:
    """A database reference of a node, with its qualifier, the same for all formats.

    Each format writes it in its own way, e.g. TabularQual as
    (qualifier, curie), SBML as an RDF CV term with the url.
    """
    namespace: str      # "bqbiol" or "bqmodel"
    qualifier: str      # e.g. "is", "isVersionOf"
    prefix: str         # identifiers.org prefix, e.g. "chebi"
    local_id: str       # e.g. "2365"

    @property
    def curie(self) -> str:
        return f"{self.prefix}:{self.local_id}"

    @property
    def url(self) -> str:
        return f"{IDENTIFIERS_ORG}/{self.curie}"


class AnnotationManager:
    """Parses <db>:<id> references into Annotations, driven by a YAML registry of databases."""

    DEFAULT_QUALIFIER = "bqbiol:isVersionOf"

    def __init__(self, yaml_filepath: Optional[Path] = None):
        self._registry: Dict[str, Dict[str, Any]] = {}
        self.load_from_yaml(yaml_filepath)

    def load_from_yaml(self, yaml_filepath: Optional[Path] = None):
        """Loads the database registry."""
        if yaml_filepath is None:
            yaml_filepath = Path(__file__).with_name("annotation_registry.yaml")

        if not yaml_filepath.exists():
            print(f"Warning: Annotation mapping config not found at '{yaml_filepath}'. Running empty.")
            return

        with open(yaml_filepath, "r", encoding="utf-8") as file:
            config = yaml.safe_load(file)

        # store the whole metadata dictionary as-is
        for prefix, meta in config.get("databases", {}).items():
            self._registry[prefix.lower()] = meta
            # ensure a canonical prefix exists
            if "canonical_prefix" not in self._registry[prefix.lower()]:
                self._registry[prefix.lower()]["canonical_prefix"] = prefix.lower()

    def process_node(self, db_references: List[str]) -> Tuple[List[Annotation], List[str]]:
        """Annotations for the references of a node, and the references that
        could not be parsed or whose database is not in the registry."""

        annotations: List[Annotation] = []
        invalid_refs: List[str] = []

        for ref in db_references:
            if ":" not in ref:
                invalid_refs.append(ref)
                continue

            prefix, local_id = ref.split(":", 1)
            meta = self._registry.get(prefix.strip().lower())

            if meta is None:
                invalid_refs.append(ref)
                continue

            namespace, qualifier = meta.get("qualifier", self.DEFAULT_QUALIFIER).split(":", 1)
            annotations.append(Annotation(
                namespace=namespace,
                qualifier=qualifier,
                prefix=meta["canonical_prefix"],
                local_id=local_id.strip(),
            ))

        return annotations, invalid_refs

# shared instance
annotation_manager = AnnotationManager()
