#!/usr/bin/env python3
'''
Definition of a taxonomy tree
'''

class TaxNode:
    '''
    Nodes belonging to a taxonomy tree
    '''
    def __init__(self, name, parent=None, children=[]):
        '''
        parent: the parent in the tree
        children: list of known children
        '''
        self.name = name
        self.parent = parent
        self.children = children

    def add_child(self, new_child):
        '''
        Add a new child to the list of children
        '''
        if new_child.name not in [child.name for child in self.children]:
            self._children.append(new_child)

    def child_end_count(self):
        '''
        Return the number of end-point nodes, i.e. the number of species 
        in the current branch of the tree.
        '''
        if len(self.children) == 0:
            return 1
        count = 0
        for child in self.children:
            count += child.child_end_count()

        return count

    @property
    def children(self):
        '''
        Return the list of children
        '''
        return self._children

    @children.setter
    def children(self, new_children):
        '''
        Change to a new list of children
        '''
        self._children = list(new_children)

    @property
    def name(self):
        '''
        Return the name of the node
        '''
        return self._name

    @name.setter
    def name(self, new_name):
        '''
        List the name of the node
        '''
        self._name = new_name

    @property
    def parent(self):
        '''
        Return the parent in the tree
        '''
        return self._parent

    @parent.setter
    def parent(self, new_parent):
        '''
        Set the parent in the tree
        '''
        self._parent = new_parent
        if type(new_parent) is type(self):
            new_parent.add_child(self)
            # if child failed to be added, eg name already existed :
            if self not in new_parent.children:
                self._parent = None
