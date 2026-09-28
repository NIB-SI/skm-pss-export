'''
Graph class for Neo4j database interaction.
This class provides methods to connect to a Neo4j database
and execute queries.

Internal to the package: PSSAdapter opens a connection for each collection
of data and closes it afterwards, callers only pass connection settings.
'''
import os
from dotenv import dotenv_values

from neo4j import GraphDatabase

# connection setting -> environment / .env variable
ENV_VARIABLES = {
    'uri': 'MY_NEO4J_URI',
    'user': 'MY_NEO4J_USER',
    'pwd': 'MY_NEO4J_PASSWORD',
}


def resolve_connection_settings(uri=None, user=None, pwd=None):
    '''
    Resolve the database connection settings. For each setting, use (in order):
        1. the argument, if not None
        2. the environment variable (MY_NEO4J_URI, MY_NEO4J_USER, MY_NEO4J_PASSWORD)
        3. the same variable in a .env file in the current folder

    Empty strings are valid values (e.g. a database without authentication).

    Returns
    -------
    dict with keys 'uri', 'user', 'pwd'
    '''
    settings = {'uri': uri, 'user': user, 'pwd': pwd}

    dotenv = dotenv_values(".env") if os.path.exists('.env') else {}

    for key, variable in ENV_VARIABLES.items():
        if settings[key] is None:
            settings[key] = os.environ.get(variable, dotenv.get(variable))

    missing = [ENV_VARIABLES[k] for k, v in settings.items() if v is None]
    if missing:
        raise ValueError(
            "Missing database connection settings, pass them as arguments or "
            f"set in the environment or a .env file: {', '.join(missing)}")

    return settings


class GraphDB:

    def __init__(self, uri, user, pwd):
        '''
        Connect to the database.
        Use as a context manager to make sure the connection is closed:

            with GraphDB(**resolve_connection_settings()) as graph_db:
                graph_db.run_query(...)
        '''

        # connect to the database
        try:
            print(f"Connecting to database at {uri} with user {user}")
            self.driver = GraphDatabase.driver(uri, auth=(user, pwd))
        except Exception as e:
            raise ConnectionError(f"Failed to connect to the database: {e}")

        # verify connection
        try:
            self.driver.verify_connectivity()
        except Exception:
            self.driver.close()
            raise
        print("Connection established.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.driver.close()

    def run_query(self, query_function, *args):

        results = []
        with self.driver.session() as session:
            query_result = session.execute_read(query_function, *args)
            # current_app.logger.info(query_result)

            for r in query_result:
                results.append(r)

        return results
