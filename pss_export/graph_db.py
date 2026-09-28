'''
Graph class for Neo4j database interaction.
This class provides methods to connect to a Neo4j database
and execute queries.
'''
import os
from dotenv import dotenv_values

from neo4j import GraphDatabase

class GraphDB:

    def __init__(self, uri=None, user=None, pwd=None):
        '''
        Initialize the Graph class.
        If no parameters are provided, look for a .env file
        load the database connection parameters from it.
        '''

        # if any of the parameters are missing, look for a .env file
        if uri is None or user is None or pwd is None:

            # if an env file exists in the current folder, load it
            if os.path.exists('.env'):
                config = dotenv_values(".env")
                if uri is None:
                    uri = config.get('MY_NEO4J_URI', None)
                if user is None:
                    user = config.get('MY_NEO4J_USER', None)
                if pwd is None:
                    pwd = config.get('MY_NEO4J_PASSWORD', None)

        # if any of the parameters are still missing, raise an error
        if uri is None or user is None or pwd is None:
            raise ValueError("Missing database connection parameters: uri, user, and pwd are required.")

        # connect to the database
        try:
            print(f"Connecting to database at {uri} with user {user}")
            self.driver = GraphDatabase.driver(uri, auth=(user, pwd))
        except Exception as e:
            raise ConnectionError(f"Failed to connect to the database: {e}")


        # verify connection
        self.driver.verify_connectivity()
        print("Connection established.")

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
