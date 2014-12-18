from argparse import ArgumentParser
import sys
import re
from ete_dev import PhyloTree, add_face_to_node, TextFace, TreeStyle, SequenceFace, random_color, SeqMotifFace
from ete_dev.treeview import layouts
from ete_dev.parser import newick
from ete_dev.ncbi_taxonomy import ncbiquery as ncbi
tracked_clades = ["Eukaryota", "Viridiplantae", "Fungi",
                     "Alveolata", "Metazoa", "Stramenopiles", "Rhodophyta",
                     "Amoebozoa", "Crypthophyta", "Bacteria",
                     "Alphaproteobacteria", "Betaproteobacteria", "Cyanobacteria",
                     "Gammaproteobacteria"]

#newick.set_float_format("%0.32g")

def spname(name):
    m = re.search('\{([^}]+)\}', name)
    if m:
        return m.groups()[0]
    else:
        return name.split('|')[0].strip().replace('_', ' ')
        taxid = name.split('.', 1)[0]

        tax2name = ncbi.get_taxid_translator([taxid])
        if int(taxid) not in tax2name:
            print 'name', name        , taxid, tax2name
        return tax2name.get(int(taxid), taxid)

def draw(t, draw_alg=True):    
    def ly(node):
        node.img_style['vt_line_width'] = 1
        node.img_style['hz_line_width'] = 1
        if node.is_leaf():
            add_face_to_node(TextFace(' (%s)' %node.name.split()[0].replace("/exon2", ""), fsize=10, fgcolor='slategrey', tight_text=False), node, 1, position='branch-right')
            add_face_to_node(TextFace(node.species, fsize=12, fgcolor='black', fstyle='italic', tight_text=False), node, 0, position='branch-right')
            c = 1
            for tname in tracked_clades:
                if tname in node.named_lineage:
                    linF = TextFace(tname, fsize=10, fgcolor='white')
                    linF.margin_left = 3
                    linF.background.color = lin2color[tname]

                    add_face_to_node(linF, node, c, position='aligned')
                    c += 1
            for n in xrange(1, 20-(c-1)):
                add_face_to_node(TextFace('', fsize=10, fgcolor='slategrey'), node, c, position='aligned')
                c+=1

            if draw_alg and 'sequence' in node.features:
                #seqFace = SequenceFace(node.sequence,"aa",13)
                seqFace = SeqMotifFace(node.sequence, [])
                # [10, 100, "[]", None, 10, "black", "rgradient:blue", "arial|8|white|domain Name"],
                motifs = []
                last_lt = None
                for c, lt in enumerate(node.sequence):
                    if lt != '-':
                        if last_lt is None:
                            last_lt = c

                        if c+1 == len(node.sequence):
                            start, end = last_lt, c
                            w = end-start
                            motifs.append([start, end, "[]", w, 13, "slategrey", "slategrey", None])
                            last_lt = None

                    elif lt == '-':
                        if last_lt is not None:
                            start, end = last_lt, c-1
                            w = end-start
                            motifs.append([start, end, "[]", w, 13, "slategrey", "slategrey", None])
                            last_lt = None

                if not motifs:
                    print node, node.sequence

                seqFace = SeqMotifFace(node.sequence, motifs,
                                       intermotif_format="line",
                                       seqtail_format="line", scale_factor=1)
                add_face_to_node(seqFace, node, 20, aligned=True)

        else:
            if node.up: 
                add_face_to_node(TextFace('% 3g' %node.support, fsize=11, fgcolor='indianred'), node, 0, position='branch-top')
                
            if hasattr(node, "support2") and node.up:
                add_face_to_node(TextFace('% 3g' %float(node.support2), fsize=11, fgcolor='steelblue'), node, 0, position='branch-bottom')

        node.img_style['size'] = 0
        node.img_style['hz_line_color'] = 'black'
        node.img_style['vt_line_color'] = 'black'

    colors = random_color(num=len(tracked_clades))
    lin2color = dict([(ln, colors[i]) for i, ln in enumerate(tracked_clades)])
    ts = TreeStyle()
    ts.draw_aligned_faces_as_table = False
    ts.draw_guiding_lines = False
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.layout_fn = ly
    print 'Rendering tree.pdf'
    t.render('tree.svg', tree_style=ts)
    t.render('tree.png', tree_style=ts)


def simplify(stop_values, stop_attr):
    def is_leaf(node):
        #if node.name in stop_clades or len(node.children) == 0:
        for v in stop_values:
            pass
        for lname in reversed(node.named_lineage):
            if lname in stop_clades:
                node.name = lname
                return True
        if len(node.children) == 0:
            return True
        
        return False

    for lf in t:
        lf.name = str(lf.taxid)
    print t.write(format=9, is_leaf_fn=is_leaf)


def reduce_trees(trees_file):
    for t in translate_ids(trees_file, outgroup_lineage="Plasmodium"):
        if not line.strip() or line.startswith('#'):
            continue
        t.dist = 0

        stop_clades = set(['Bacteria', 'Eukaryota', "Viridiplantae", "Coccidia", "Apicomplexa", "Plasmodium"])
        
        def is_leaf(node):
            #if node.name in stop_clades or len(node.children) == 0:
            for lname in reversed(node.named_lineage):
                if lname in stop_clades:
                    node.name = lname
                    return True
            if len(node.children) == 0:
                return True
            return False

        for lf in t:
            lf.name = str(lf.taxid)
        t.unroot()
        t.sort_descendants()
        
        print t.write(format=9, is_leaf_fn=is_leaf)

def translate_ids(trees_file, outgroup_lineage="Bacteria"):
    for line in open(trees_file):
        if not line.strip() or line.startswith('#'):
            continue

        t = PhyloTree(line, sp_naming_function=spname)
        #t.set_outgroup(t.get_midpoint_outgroup())

        for lf in t:
            lf.add_features(coded_name = lf.name)            
            if lf.name in NAME2SP:
                lf.name = "%s {%s}" %(lf.name, NAME2SP[lf.name])
         
        t.dist = 0
        ncbi.connect_database()
        name2sp = ncbi.get_name_translator(t.get_species())
        for lf in t.iter_leaves():
            lf.add_features(taxid=name2sp.get(lf.species, 0))

        t.set_outgroup(t.search_nodes(taxid=9606)[0])
        ncbi.annotate_tree(t, attr_name='taxid')
        t.set_outgroup(t.get_common_ancestor([lf for lf in t if outgroup_lineage in lf.named_lineage]))
        ncbi.annotate_tree(t, attr_name='taxid')
            
        #print t.write(features=[])
        #print t.write()
        yield t


def dump_constraints(t):
    leaves = set(t.get_leaves())
  
    rho = set([lf for lf in leaves if "Rhodophyta" in lf.named_lineage])
    stra = set([lf for lf in leaves if "Stramenopiles" in lf.named_lineage])
    plants = set([lf for lf in leaves if "Viridiplantae" in lf.named_lineage])
    
    alv = set([lf for lf in leaves if "Alveolata" in lf.named_lineage])
    euk = set([lf for lf in leaves if "Eukaryota" in lf.named_lineage]) 
    meta = set([lf for lf in leaves if "Metazoa" in lf.named_lineage])
    

    plas = set([lf for lf in leaves if "Plasmodium" in lf.named_lineage])
    api = set([lf for lf in leaves if "Apicomplexa" in lf.named_lineage])
    bact = set([lf for lf in leaves if "Bacteria" in lf.named_lineage])

    c = 0
    for tname, target in [("Api", api)]:
        for gname, group in [("Rho", rho), ("Stra", stra), ("Alg", stra|rho), ("Alv", alv)]:#, ("Plants", plants)]:
            group = group - target 
            for size in xrange(1, len(group)+1):
                for cname, comb in enumerate(it.combinations(group, size)):
                    c += 1
                    hyp = set(comb)|target
                    rest = leaves - hyp
                    OUTFILE = "constraint_%s+%s_VS_rest" %(tname, "%s-%s-%s" %(gname, size, cname))
                    open("constraints/"+OUTFILE, 'w').write("(%s, ((%s),%s));\n" %(
                            ','.join(map(lambda lf: lf.coded_name, rest)),
                            ','.join(map(lambda lf: lf.coded_name, target)),
                            ','.join(map(lambda lf: lf.coded_name, set(comb)))))
                    cmd = "/home/huerta/standard-RAxML-master/raxmlHPC-SSE3 -m PROTGAMMAWAG -f d -g /g/bork/huerta/prx5/prx5/constraints/%s -n %s -s /g/bork/huerta/prx5/prx5/manual.fa -w /g/bork/huerta/prx5/prx5/raxml_constraint_test/ -p 31416;" %(OUTFILE, OUTFILE)
                    print cmd
    print >>sys.stderr, c 

def content_id(node):
    import hashlib
    return hashlib.md5(str(sorted(node.get_leaf_names()))).hexdigest()
    
    
        
if __name__ == '__main__':
    from string import strip
    import itertools as it
    NAME2SP = dict([map(strip, line.split('\t')) for line in open('gname2sp')])

    parser = ArgumentParser()
    parser.add_argument('--reduce', dest='reduce', type=str)
    parser.add_argument('--translate', dest='translate', type=str)
    parser.add_argument('--constraints', dest='constraints', type=str)
    parser.add_argument('--draw', dest='draw', type=str)
    parser.add_argument('--combine', dest='combine', type=str, nargs="+")
    args = parser.parse_args()


    
    if args.combine:
        target = translate_ids(args.combine[0]).next()
        source = translate_ids(args.combine[1]).next()
        print set(source.get_leaf_names())  ^ set(target.get_leaf_names())
        tid2node = {}
        for n in target.traverse():
            tid2node[content_id(n)] = n
            
        for n in source.traverse():
            sid = content_id(n)
            if sid in tid2node:
                tid2node[sid].add_features(support2 = n.support)

        print target.write(features=['support2'])
        target.ladderize()
        draw(target, False)
    
    if args.translate:
        for t in translate_ids(args.translate):
            for lf in t:
                lf.name = lf.taxid
                lf.name = '_'.join(lf.spname.split()[:2]).replace('.', "").replace("/", "")
            print t
            print t.write()

    if args.reduce:
        reduce_trees(args.reduce)

    if args.constraints:
        for t in translate_ids(args.constraints):
            dump_constraints(t)

    if args.draw:
        for t in translate_ids(args.draw):
            t.ladderize()
            draw(t, False)


