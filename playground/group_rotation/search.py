#!/usr/bin/python

import collections

def find_children(root, blocked):
    """
    Finds all the children of a given atom on one side of the bond.
    """
    # Create a queue for the fringe (i.e. breadth-first search).
    fringe = collections.deque()
    # Add root to the fringe.
    fringe.appendleft(root)
    # Add root onto the visited list
    visited = [root]
    # While there are things on the fringe, pop (from the right) of
    # the fringe, examine neighbouring nodes and then add them to the
    # fringe if they haven't been visited or blocked.
    while len(fringe) > 0:
        current = fringe.pop()
        for node in current.bonded:
            if (node not in visited) and (node not in blocked):
                visited.append(node)
                fringe.appendleft(node)
    return visited
