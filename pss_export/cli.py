#!/usr/bin/env python3

import click
import functools

from .pss.config import pss_schema_config

SPECIES = list(pss_schema_config.species)

from pss_export import PSSAdapter

# click option that converts comma separated string into list
# if no argument is provided, it returns None (instead of empty list)
# (existing options in click do not handle unlimited number of values)
# simple solution found here: https://stackoverflow.com/a/48394085
class ConvertStrToList(click.Option):
    def type_cast_value(self, ctx, value):
        try:
            if value is None:
                return None
            else:
                return [v.strip() for v in value.split(",")]
        except Exception:
            raise click.BadParameter(value)


def neo4j_common_params(func):
    @click.option("--neo4j-uri", default=None, help="Neo4j connection URI (default: MY_NEO4J_URI from the environment or .env).")
    @click.option("--neo4j-user", default=None, help="Neo4j username (default: MY_NEO4J_USER from the environment or .env).")
    @click.option("--neo4j-password", default=None, help="Neo4j password (default: MY_NEO4J_PASSWORD from the environment or .env).")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper

def export_common_params(func):
    # model_id, model_name, model_description, etc
    @click.option("--model-id", default="my_pss_model", help="Model ID to use in export.")
    @click.option("--model-name", default="PSS Model", help="Model name to use in export.")
    @click.option("--model-description", default="Model exported from the Plant Stress Signalling knowledge graph (PSS) available at https://skm.nib.si using the skm-pss-export package.", help="Model description to use in export.")
    @click.option("--model-version", default=None, help="Version of the exported model, e.g. the PSS version.")
    @click.option("--creator", default=None, help="Creator of the model, in the format of: familyName | givenName | organization | email. Can be specified multiple times for multiple creators.", multiple=True)
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper

def modelfixing_common_params(func):
    @click.option("--model-fixes-identify", is_flag=True, help="Run the model fixing module to identify model inconsistencies and suggest automatic fixes.")
    @click.option("--model-fixes-apply", is_flag=True, help="If running the model fixing module, also apply all suggested model fixes")
    @click.option("--model-fixes-interactive", is_flag=True, help="If running the model fixing module, enter interactive mode to visualise and optionally apply fixes per node. Overrides `--model-fixes-apply`. ")
    @click.option("--nodes-to-ignore", default='default', help="Node to leave out of the model (SBML, TabularQual). Default is 'default': the nodes_to_ignore in pss_export_config.yaml.")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper

def reaction_filter_common_params(func):
    @click.option("--reactions", cls=ConvertStrToList, default=None, help="Comma-separated list of reaction IDs to include in export.")
    @click.option("--access", default='public', type=click.Choice(['public', 'restricted', 'all']), help="Data access level: 'public' (default) or 'restricted' (all reactions; 'all' is an alias).")
    @click.option("--species", default='ath', type=click.Choice(SPECIES + ['all']), callback=lambda ctx, param, value: None if value == 'all' else value,
                  help="Limit to the reactions whose functional clusters all have genes in this species (default: ath); gene annotations (tair: for ath) and the gene network are for this species. 'all': no species filter (not for the gene network).")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper

@click.group()
def cli():
    """Top-level CLI entrypoint."""
    pass

@cli.command()
@neo4j_common_params
@export_common_params
@reaction_filter_common_params
@modelfixing_common_params
@click.argument("filename", type=click.Path())
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
@click.option("--include-genes", is_flag=True, help="Include genes as species (transcription reactions).")
@click.option("--entities-table", default=None, type=click.Path(), help="Path to also export a table of entities in model.")
@click.option("--kinetic-laws", is_flag=True, help="Include kinetic laws (SBO term only) in SBML output.")
def to_sbml(neo4j_uri, neo4j_user, neo4j_password,
            model_id, model_name, model_description, model_version, creator,
            access, species, reactions,
            model_fixes_identify, model_fixes_apply, model_fixes_interactive,
            nodes_to_ignore,
            filename,
            verbose,
            include_genes,
            entities_table,
            kinetic_laws):
    """
    Export model to SBML.
    """
    if verbose:
        click.echo(f"Exporting to SBML...")
        click.echo(f"  Access: {access}")
        click.echo(f"  Neo4j URI: {neo4j_uri}")
        click.echo(f"  Output file: {filename}")

    try:

        # build adapter (connects to PSS in neo4j only while collecting)
        adapter = PSSAdapter(neo4j_uri=neo4j_uri,
                    neo4j_user=neo4j_user,
                    neo4j_password=neo4j_password,
                    model_id=model_id,
                    model_name=model_name,
                    model_description=model_description,
                    model_version=model_version,
                    creator=creator)
        adapter.collect_reactions(reactions=reactions, access=access, species=species, include_genes=include_genes, nodes_to_ignore=nodes_to_ignore)
        if model_fixes_identify:
            adapter.model_fixes(apply_fixes=model_fixes_apply, interactive=model_fixes_interactive)

        adapter.create_sbml(filename=filename,
                    entities_table=entities_table,
                    kinetic_laws=kinetic_laws)

        if verbose:
            click.secho("SBML export complete.", fg="green")

    except Exception as e:
        click.secho(f"Error: {e}", fg="red", err=True)
        raise click.Abort()


@cli.command()
@neo4j_common_params
@export_common_params
@reaction_filter_common_params
@modelfixing_common_params
@click.argument("filename", type=click.Path())
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def to_tabularqual(neo4j_uri, neo4j_user, neo4j_password,
            model_id, model_name, model_description, model_version, creator,
            access, species, reactions,
            model_fixes_identify, model_fixes_apply, model_fixes_interactive,
            nodes_to_ignore,
            filename,
            verbose
            ):
    """
    Export model to TabularQual.
    """

    if verbose:
        click.echo(f"Exporting to TabularQual...")
        click.echo(f"  Access: {access}")
        click.echo(f"  Neo4j URI: {neo4j_uri}")
        click.echo(f"  Output file: {filename}")

    # build adapter (connects to PSS in neo4j only while collecting)
    adapter = PSSAdapter(neo4j_uri=neo4j_uri,
                    neo4j_user=neo4j_user,
                    neo4j_password=neo4j_password,
                    model_id=model_id,
                    model_name=model_name,
                    model_description=model_description,
                    model_version=model_version,
                    creator=creator)
    adapter.collect_reactions(reactions=reactions, access=access, species=species, nodes_to_ignore=nodes_to_ignore)
    if model_fixes_identify:
        adapter.model_fixes(apply_fixes=model_fixes_apply, interactive=model_fixes_interactive)

    adapter.create_tabularqual(filename=filename)

    click.echo(f"Wrote spreadsheet to {filename}")

@cli.command()
@neo4j_common_params
@export_common_params
@reaction_filter_common_params
@modelfixing_common_params
@click.argument("filename", type=click.Path())
@click.option("--newt", is_flag=True, help="The Newt variant: with Newt's custom properties and colours.")
def to_sbgn(neo4j_uri, neo4j_user, neo4j_password,
            model_id, model_name, model_description, model_version, creator,
            access, species, reactions,
            model_fixes_identify, model_fixes_apply, model_fixes_interactive,
            nodes_to_ignore,
            filename, newt):
    """
    Export the model reactions to SBGN-ML (Process Description), laid out with Graphviz.
    """
    adapter = PSSAdapter(neo4j_uri=neo4j_uri, neo4j_user=neo4j_user, neo4j_password=neo4j_password,
                         model_id=model_id, model_name=model_name, model_description=model_description,
                         model_version=model_version, creator=creator)
    adapter.collect_reactions(reactions=reactions, access=access, species=species, nodes_to_ignore=nodes_to_ignore)
    if model_fixes_identify:
        adapter.model_fixes(apply_fixes=model_fixes_apply, interactive=model_fixes_interactive)
    n = adapter.create_sbgn_newt(filename) if newt else adapter.create_sbgn(filename)
    click.echo(f"Wrote {n} reactions to {filename}")


def _collect_for_network(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions):
    # the networks keep every node: nodes_to_ignore is for the models only
    adapter = PSSAdapter(neo4j_uri=neo4j_uri, neo4j_user=neo4j_user, neo4j_password=neo4j_password)
    adapter.collect_reactions(reactions=reactions, access=access, species=species, nodes_to_ignore=None)
    return adapter


def network_params(func):
    @neo4j_common_params
    @reaction_filter_common_params
    @click.argument("edges_file", type=click.Path())
    @click.argument("nodes_file", type=click.Path())
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper


@cli.command()
@network_params
def to_reaction_graph(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions, edges_file, nodes_file):
    """
    Export the reaction graph (extended SIF edges and nodes): entities and reactions.
    """
    adapter = _collect_for_network(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions)
    n = adapter.create_reaction_graph(edges_file, nodes_file)
    click.echo(f"Wrote {n} edges to {edges_file}, nodes to {nodes_file}")


@cli.command()
@neo4j_common_params
@reaction_filter_common_params
@click.argument("filename", type=click.Path())
def to_reaction_graph_json(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions, filename):
    """
    Export the reaction graph as one JSON file (nodes and edges): the data of the PSS Explorer
    (use --species all).
    """
    adapter = _collect_for_network(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions)
    n = adapter.create_reaction_graph_json(filename)
    click.echo(f"Wrote {n} edges to {filename}")


@cli.command()
@network_params
def to_interaction_network(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions, edges_file, nodes_file):
    """
    Export the interaction network (extended SIF edges and nodes): entity -> entity influences.
    """
    adapter = _collect_for_network(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions)
    n = adapter.create_interaction_network(edges_file, nodes_file)
    click.echo(f"Wrote {n} edges to {edges_file}, nodes to {nodes_file}")


@cli.command()
@network_params
def to_gene_network(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions, edges_file, nodes_file):
    """
    Export the gene network (extended SIF edges and nodes): the interaction network with functional
    clusters expanded into their genes in --species (not 'all').
    """
    if species is None:
        raise click.UsageError("A gene network needs a species (--species), not 'all'.")
    adapter = _collect_for_network(neo4j_uri, neo4j_user, neo4j_password, access, species, reactions)
    n = adapter.create_gene_network(edges_file, nodes_file)
    click.echo(f"Wrote {n} edges to {edges_file}, nodes to {nodes_file}")


@cli.command()
@neo4j_common_params
@click.option("--access", default="public", type=click.Choice(["public", "restricted"]), help="Access level of the reactions described.")
@click.argument("filename", type=click.Path())
def to_faidare(neo4j_uri, neo4j_user, neo4j_password, access, filename):
    """
    Export the FAIDARE data discovery file (JSON): an entry per gene and species.
    """
    adapter = PSSAdapter(neo4j_uri=neo4j_uri, neo4j_user=neo4j_user, neo4j_password=neo4j_password)
    adapter.collect_reactions(access=access, species=None, nodes_to_ignore=None)
    n = adapter.create_faidare(filename)
    click.echo(f"Wrote {n} entries to {filename}")


@cli.command()
def formats():
    """
    List the export formats.
    """
    from .formats import FORMATS
    for fmt in FORMATS.values():
        files = ", ".join(f"{f.name} (.{f.extension})" for f in fmt.files)
        click.echo(f"{fmt.key}: {fmt.title} -- {files}")


if __name__ == "__main__":
    cli()
