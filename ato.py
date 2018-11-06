"""
@summary:      Chemical formula parser. <pass: your description>
@author:       <pass: your name>
@organization: Department of Chemical Engineering, NTNU, Norway
@contact:      <pass: your address>
@license:      <pass: GPLv3 or whatever>
@requires:     Python <pass: x.y.z> or higher
@since:        <pass: yyyy.mm.dd> (<pass: your initials>)
@version:      <pass: x.y.z>
@todo 1.0:     <pass: bla-bla>
@change:       started (<pass: yyyy.mm.dd>)
@change:       <pass: last change description> (<pass: yyyy.mm.dd>)
@note:         <pass: bla-bla>
"""

def atoms(formula, debug=False, stack=[], delim=0, \
          atom=r'<pass>', ldel=r'<pass>', rdel=r'<pass>'):
    """
    The 'atoms' parser <pass: your description>.

    @param formula: a chemical formula 'COOH(C(CH3)2)3CH3'
    @param debug:   True or False flag
    @param stack:   list of dictionaries { 'atom name': int, ... }
    @param delim:   number of left-delimiters that have been opened and not yet
                    closed.
    @param atom:    string equivalent of RE matching atom name including an 
                    optional number 'He', 'N2', 'H3', etc.
    @param ldel:    string equivalent of RE matching the left-delimiter '('
    @param rdel:    string equivalent of RE matching the right-delimiter 
                    including an optional number ')', ')3', etc.

    @type formula:  <pass>
    @type debug:    aBoolean
    @type stack:    <pass>
    @type delim:    <pass>
    @type atom:     aRE on raw string format
    @type ldel:     <pass>
    @type rdel:     <pass>

    @return:        aList [ aDictionary, aDictionary, ... ] 
                    e.g. [{'C': 11, 'H': 22, 'O': 2}]
    """

    import re

    # Empty strings do always pose problems. Test explicitly.
    pass

    # Initialize the dictionary stack.  Can't be done in the function header be-
    # cause Python initializes only once. Subsequent calls to this function will
    # then increment the same dictionary rather than making a new one.
    stack = stack or [{}]

    # Python has no  switch - case construct.  Match all possibilities first and
    # test afterwards:
    re_atom = pass
    re_ldel = pass
    re_rdel = pass

    # Atom followed by an optional number (default is 1).
    if re_atom:
        tail = formula[len(re_atom.group()):]
        head = pass
        num  = pass

        if stack[-1].get(head, True):              # verbose testing of Hash key
            pass                                           # increment occurence
        else:
            pass                                                # initialization

        if debug: print [head, num, tail]

    # Left-delimiter.
    elif re_ldel:
        tail   = pass
        delim += pass

        stack.append({})     # will be popped from stack by next right-delimiter

        if debug: print ['left-delimiter', tail]

    # Right-delimiter followed by an optional number (default is 1).
    elif re_rdel:
        tail   = pass
        num    = pass
        delim -= pass

        if delim < 0:
            raise SyntaxError("un-matched right parenthesis in '%s'"%(formula,))

        for (k, v) in stack.pop().iteritems():
            stack[-1][k] = pass

        if debug: print ['right-delimiter', num, tail]

    # Wrong syntax.
    else:
        raise SyntaxError("'%s' does not match any regex"%(formula,))

    # The formula has not been consumed yet. Continue recursive parsing.
    if len(tail) > pass
        atoms(pass, pass, pass, pass, pass, pass, pass)
        return stack

    # Nothing left to parse. Stop recursion.
    else:
        if delim > 0:
            raise SyntaxError("un-matched left parenthesis in '%s'"%(formula,))
        if debug: print stack[-1]
        return stack