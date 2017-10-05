# InventorAnalyze
# To construct a query using the PatentsView API, see: www.patentsview.org/api/query-language.html
# Be mindful of the default number of results returned by PatentsView

from argparse import ArgumentParser
from collections import defaultdict
from json import load as jsonload
from math import ceil
from textwrap import wrap
from urllib import urlopen
from warnings import catch_warnings, simplefilter

from matplotlib import pyplot
from networkx import Graph, DiGraph, degree_assortativity_coefficient, bipartite, eigenvector_centrality, pagerank, \
    gnm_random_graph, articulation_points, draw_networkx_edges, connected_components, density, draw_networkx, \
    shell_layout, connected_component_subgraphs
from networkx.utils.rcm import cuthill_mckee_ordering


_help_message = """
To analyze a network of patents and inventors, pass a URL as an argument using the PatentsView API:
See http://www.patentsview.org/api/patent.html for the API query structure.
""".strip()

_further_help = """
For example, enter the following using single quotes for the URL argument:
python InventorAnalyze.py 'http://www.patentsview.org/api/patents/query?q={"cpc_subgroup_id":"C12Q1/6886"}&o={"per_page":5000}'
""".strip()

_assortativity_explanation = """The disambiguation of the names of joint inventors yields a social 
network of network of inventors wherein the edges between the nodes represent a joint inventorship. 
Some of these inventors have a high degree of collaboration because they have 
invented with a myriad of co-inventors since 1976; others have co-invented with
only a few other inventor since 1976.  How correlated are the collaboration degrees
of the inventors in these co-inventor pairs?  For example, do high collaborators work with 
high collaborators and low with low?  We measure this correlation using the graph-level 
property called degree assortativity, which ranges from -1 to +1.  The assortativity measure
is positive when there is a positive correlation between the collaboration degrees across the graph.""".strip()

_pagerank_explanation = """Patents:  Top 10 pageranks based on citations (calculated by disambiguating inventors 
to remove their self-citations.)  To see all of the pageranks, run again using the -p parameter.""".strip()

_eigenvector_centrality_explanation = """Inventors:  Top 10 eigenvector centralities, calculated by first disambiguating inventors 
and then weighing each of their collaborations by the total of the pageranks of their co-invented 
patents.  To see all of the centralities, run again using the -e parameter.""".strip()

def main():
    def the_weight(GP, u, v):
        '''
        This will be used to weights of the co-inventor graph in a way that is weighted by the pagerank.
        :param GP: the graph from which the weights are being computed
        :param u: a first node
        :param v: a second node
        :return: the weight for the edge between the nodes u and v
        '''
        uset = set(GP[u])
        vset = set(GP[v])
        x = sum([pr.get(patent, 0.0) for patent in
                 uset.intersection(vset)])  # pr will be computed before this function is called
        return x

    def make_str(input_string):
        '''
        
        :param input_string: a string that potentially has unprintable characters
        :return: a string omitting unprintable characters
        '''
        result = "".join([letter for letter in input_string if ord(letter) < 128])
        return result

    def output_assort(graphassort):
        '''
        
        :param graphassort: a graph for which assortativity should be calculated
        :return: a printable string indicating the assortativity of the graph and the
                    significance of the assortativity compared to random graphs using 
                    the conditional uniform graph test.
        '''
        if density(graphassort) == 1.0:
            result = "Complete Graph"
        else:
            nbootstrap = 200
            this_assortativity = degree_assortativity_coefficient(graphassort)
            random_distributions_this_graph = [degree_assortativity_coefficient(
                gnm_random_graph(graphassort.number_of_nodes(), graphassort.number_of_edges())) for _ in
                                               xrange(nbootstrap)]
            this_p = float(
                len([datum for datum in random_distributions_this_graph if datum > this_assortativity])) / float(
                len(random_distributions_this_graph))
            result = "\nDegree Assortativity: " + str(float(round(1000 * this_assortativity)) / 1000.0) + " (" + str(
                100.0 - float(round((this_p * 100000.0) / 1000.0))) + " percentile)"
        return result

    show_graphs = True  # show visualizations 1 and 2
    n = 16  # number of graphs to display in a grid in the first visualization
    top_patent_number = 250  # number of edges in the second visualization

    help_message = _help_message
    further_help = _further_help
    parser = ArgumentParser(description=help_message)
    parser.add_argument('url', metavar='URL', type=str, nargs='*', help=further_help)
    parser.add_argument("-e", "--all_eigenvector_centralities", action="store_true")
    parser.add_argument("-p", "--all_pageranks", action="store_true")
    args = parser.parse_args()

    if args.all_eigenvector_centralities:
        output = "e"
    elif args.all_pageranks:
        output = "p"
    else:
        output = "a"

    if args.url == []:
        print help_message
        print further_help
    else:
        query = args.url[
                    0] + '&f=["patent_number","inventor_id","inventor_first_name","inventor_last_name",' \
                         '"cited_patent_number","patent_title"]'
        # retrieve these fields in addition to what was specified by the user's URL

        BP = Graph()  # A BiPartite Graph (BP) of applications and inventors
        DG = DiGraph()  # A Directed Graph for the citation network
        inventors = defaultdict(frozenset)
        # a dictionary of applications mapping to their inventors to eliminate self-citations
        titles = defaultdict(str)
        try:
            data = jsonload(urlopen(query))
            assert not (data[u'patents'] is None)
            queryset = set([])  # for later filtering the network to include only queried patents
            allinventors = set([])  # for later computing the co-inventor network
            for patent in data[u'patents']:
                patent_number = patent[u'patent_number']
                if patent_number[0] != "D":  # exclude design patents
                    queryset.add(patent_number)  # for later filtering the network to include only queried patents
                    inventorset = set(
                        [])  # for keeping track of an application's inventors to later remove self-citations
                    for inventor in patent[u'inventors']:
                        inventor_info = (
                        inventor[u'inventor_id'], inventor[u'inventor_first_name'], inventor[u'inventor_last_name'])
                        BP.add_edge(patent_number,
                                    inventor_info)
                        # A BiPartite Graph (BP) with patents and inventors, later to be projected into a co-inventor network
                        inventorset.add(
                            inventor_info)
                        # for keeping track of an application's inventors to remove self-citations later
                    inventorset = frozenset(inventorset)
                    inventors[
                        patent_number] = inventorset  # for later when removing self-citations from the citation network
                    allinventors.update(inventorset)
                    titles[patent_number] = patent[u'patent_title']
                    for citation in patent[u'cited_patents']:  # for constructing the citation network
                        citation_number = citation[u'cited_patent_number']
                        if citation_number is not None:
                            if citation[u'cited_patent_number'][0] != "D":  # exclude design patents
                                DG.add_edge(patent_number, citation[
                                    u'cited_patent_number'])
                                # some of these citations are not in the original application query; these will later be deleted
            for patent, cited_patent in DG.edges():
                # to eliminate citations that were not in the original query and to eliminate self-citations
                ok = (patent in queryset) and (cited_patent in queryset) and (
                len(inventors[patent].intersection(inventors[cited_patent])) == 0)
                if not ok:
                    DG.remove_edge(patent, cited_patent)
            pr = pagerank(DG)
            proutput = sorted(pr.items(), key=lambda x: x[1], reverse=True)
            CoInventors = bipartite.generic_weighted_projected_graph(BP, allinventors, weight_function=the_weight)
            try:
                e = sorted(eigenvector_centrality(CoInventors, weight='weight', max_iter=10000).items(),
                           key=lambda x: x[1], reverse=True)
                eflag = True
            except:
                eflag = False
                pass
            if output == 'a':
                if len(CoInventors) > 1:
                    print
                    print str(CoInventors.number_of_nodes()) + " disambiguated inventors."
                    print str(len(data[u'patents'])) + " patents."
                    print
                    actual_assortativity = degree_assortativity_coefficient(CoInventors)
                    print _assortativity_explanation
                    print
                    print "The Degree Assortativity of this Subnetwork of Disambiguated Inventors:  " + str(
                        actual_assortativity)
                    print

                    print _pagerank_explanation
                    print "Patent Number, Pagerank"
                    for patent, pagerank_value in proutput[0:10]:
                        print "US" + str(patent) + "," + str(
                            pagerank_value) + ", http://www.patentsview.org/web/#detail/patent/" + str(
                            patent.encode('utf-8'))
                    print
                    if eflag:
                        print _eigenvector_centrality_explanation
                        print "Inventor , Centrality, Inventor Link"
                        for inventor, eigenvector_centrality_value in e[0:10]:
                            print inventor[1].encode('utf-8') + " " + inventor[2].encode(
                                'utf-8') + " , " + str(
                                eigenvector_centrality_value) + " , http://www.patentsview.org/web/#detail/inventor/" + str(
                                inventor[0])
                        print
                    else:
                        print "Eigenvector Centrality could not be computed for this graph."
                        print

                    # Visualization One
                    some_colors = ['g', 'y', 'c', 'm', 'b', 'g', 'r']
                    ord_overall = list(cuthill_mckee_ordering(CoInventors))  # ordering to use to arrange graphs
                    if show_graphs:
                        rows = ceil(n ** 0.5)
                        edict = dict(e)
                        largest_subgraphs = sorted([a_graph for a_graph in connected_component_subgraphs(CoInventors) if
                                                    (density(a_graph) < 1.0 and a_graph.number_of_nodes() > 3)],
                                                   key=lambda x: len(x), reverse=True)[0:n]
                        if len(largest_subgraphs) > 0:
                            subgraphs = [a_subgraph for a_subgraph in
                                         sorted(largest_subgraphs, key=lambda x: degree_assortativity_coefficient(x),
                                                reverse=True)][0:n]
                            for subpl, this_graph in enumerate(subgraphs):
                                this_ord = [list(cuthill_mckee_ordering(this_graph))]
                                e_dict = dict([(node, ec) for node, ec in e if node in set(this_ord[0])])
                                top_pos = shell_layout(this_graph, nlist=this_ord)
                                pyplot.figure(num=n, figsize=(40, 15), facecolor='w')
                                sp = pyplot.subplot(rows, rows, subpl + 1)
                                the_color = some_colors[subpl % len(some_colors)]
                                try:
                                    top_inventors = sorted(this_graph.nodes(), key=lambda x: edict[x], reverse=True)
                                    top_inventor = top_inventors[0]
                                    flatten_me = [BP[inv].keys() for inv in top_inventors if len(BP[inv].keys()) > 0]
                                    these_pageranks = [pr.get(pat, 0) for pat in
                                                       set([item for sublist in flatten_me for item in sublist])]
                                    top_group_patent = \
                                    sorted(list(set([item for sublist in flatten_me for item in sublist])),
                                           key=lambda x: pr.get(x, 0), reverse=True)[0]
                                except:
                                    raise
                                    top_inventor = ["", "", ""]
                                    top_group_patent = ""
                                max_length = 48
                                a_caption = "\n".join(
                                    wrap(titles[top_group_patent], max_length)[0:2]) + "...\ninvolved " + str(
                                    len(set(top_inventors))) + " collaborators (e.g., " + top_inventor[1].encode(
                                    'utf-8') + " " + top_inventor[2].encode(
                                    'utf-8') + ")\n" + "Highest Pagerank: " + str(max(these_pageranks)) + output_assort(
                                    this_graph)
                                sp.set_title(a_caption, color=the_color)
                                pyplot.axis('off')
                                pyplot.tight_layout()
                                threshold = .50
                                with catch_warnings():
                                    simplefilter("ignore")
                                    draw_networkx(this_graph, node_size=5, nodelist=e_dict.keys(), node_color=the_color,
                                                  pos=top_pos, with_labels=False, edge_color=the_color)
                            with catch_warnings():
                                simplefilter("ignore")
                                pyplot.show()

                            # Visualization Two
                            P = Graph()
                            for node1, node2, weightdict in sorted(largest_subgraphs[0].edges(data=True),
                                                                   key=lambda x: x[2]["weight"], reverse=True)[
                                                            :top_patent_number]:
                                P.add_edge(node1, node2)

                            ord1 = list(cuthill_mckee_ordering(P))
                            if ord1[0] > ord1[-1]:
                                ord1 = list(reversed(ord1))
                            blank_nodes = ["blank_node"]
                            visualized_graph = P.copy()
                            visualized_graph.add_nodes_from(blank_nodes)
                            the_pos = shell_layout(visualized_graph, nlist=[ord1 + blank_nodes])
                            a = list(set(articulation_points(P)))
                            non_a = (set(ord1) - set(a))
                            the_labels = {
                            (invid, first_name, last_name): (make_str(first_name) + " " + make_str(last_name)) for
                            invid, first_name, last_name in P.nodes()}
                            articulation_labels = {nodeid: "\n".join(name.split(" ")) for nodeid, name in
                                                   the_labels.items() if (nodeid in a)}
                            non_articulation_labels = {nodeid: "".join([partname[0] for partname in name.split(" ")])
                                                       for nodeid, name in the_labels.items() if not (nodeid in a)}
                            pyplot.axis('off')
                            pyplot.tight_layout()
                            draw_networkx_edges(P, pos=the_pos, alpha=.8)
                            draw_networkx(P, edgelist=[], nodelist=non_a, node_color='g', node_size=300, font_size=8,
                                          pos=the_pos, labels=non_articulation_labels, alpha=.8)
                            draw_networkx(P, edgelist=[], nodelist=a, node_color='c', node_size=2000, font_size=8,
                                          pos=the_pos, labels=articulation_labels, alpha=.8)
                            pyplot.show()
                        else:
                            print "The network is too small to visualize properly."
                else:
                    print "The network is too small."
            elif output == 'e':
                for inventor, eigenvector_centrality_value in e:
                    print inventor[1].encode('utf-8') + "," + inventor[2].encode('utf-8') + "," + str(
                        eigenvector_centrality_value) + ", http://www.patentsview.org/web/#detail/inventor/" + str(
                        inventor[0])
            elif output == 'p':
                for patent, pagerank_value in proutput:
                    print patent.encode('utf-8') + "," + str(
                        pagerank_value) + ", http://www.patentsview.org/web/#detail/patent/" + str(
                        patent.encode('utf-8'))
        except (IOError, TypeError, ValueError):
            print "The URL that you entered probably has a typo."
            raise
        except AssertionError:
            print "The URL that you entered did not retrieve results."


if __name__ == "__main__":
    main()
