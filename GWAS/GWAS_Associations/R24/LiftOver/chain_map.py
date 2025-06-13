import re
from typing import Any, Union
from collections import defaultdict

######### Exceptions ##########
class NodeNotFound(Exception):
    pass

class PathNotFound(Exception):
    pass

### Functions ###
def find_shortest_path(map, start: str, end: str) -> list[Any]:
    dist = {start: [start]}
    queue = [start]
    
    while queue:
        current = queue.pop()
        
        if current not in map: continue

        for next in map[current]:
            if next not in dist:
                dist[next] = dist[current] + [next]
                queue.append(next)

    return dist.get(end, [])


################ Chain Files ###################
class ChainFiles:
    def __init__(self, files: dict[Any, dict[Any, str]]):
        self.files = files

    def get(self, source: str, target: str) -> str:
        return self.files[source][target]

    def graph(self):
        return {k: [k2 for k2 in self.files[k].keys()] for k in self.files.keys()}

    def find_path(self, start: str, end: str):
        if path := find_shortest_path(self.graph(), start, end):
            return path
        raise PathNotFound(f"Path not found: {start} -> {end}")

    @property
    def nodes(self):
        return list(self.graph().keys())

    def find_all_paths_to(self, target: str):
        if not target in self.nodes:
            raise NodeNotFound(f"Node {target} not found in graph...")

        candidates = self.nodes
        candidates.remove(target)

        found = []
        for candidate in candidates:
            try:
                if self.find_path(candidate, target):
                    found.append(candidate)
            except PathNotFound:
                continue
        
        return found

    @classmethod
    def from_names(cls, names: Union[str, list[str]]):
        chain_files = defaultdict(dict)

        if isinstance(names, str):
            names = [names]

        for file in names:
            match = re.findall("hg\d\d", file)

            if len(match) != 2:
                raise Exception("Could not parse file name...")
            
            chain_files[match[0]][match[1]] = file

        return cls(chain_files)