import spherogram as sph
from spherogram.links import Link, Strand, Crossing

def flype(link : Link, in2 : int, out1 : int, out2 : int):
    """
    Returns the flyped version of a link, where the tangle to be flyped is specified by incoming strands in1 = 0, in2 and outgoing strands out1, out2. 

    Parameters:
        link: the link to be flyped. 
            - The link must have 1 component (i.e. be a knot)
            - The PD presentation must correspond to traversing the link starting from the edge on the under-strand just before the crossing adjacent to the tangle and moving towards the tangle, with that starting edge labeled 0. 
            - The strands in the crossing adjacent to the tangle both must orient towards the tangle.

        in2: the label of the edge on the over-strand just before the crossing adjacent to the tangle.

        out1: the label of the edge leaving the tangle given that the tangle was entered from the under-strand of the adjacent crossing.

        out2: the label of the edge leaving the tangle given that the tangle was entered from the over-strand of the adjacent crossing.
    """
    original_pd = link.PD_code()
    num_edges = 1 + max([edge for tuple in original_pd for edge in tuple])
    flyped_pd = []
    # For-loop deals with all crossings that are not the initial one (which is flipped to the other side of the tangle)
    for i in range(1, len(original_pd)):
        # If crossing outside of tangle, doesn't change
        if (original_pd[i][0] >= out1 and original_pd[i][0] < in2) or (out2 != 0 and original_pd[i][0] >= out2):
            flyped_pd.append(list(original_pd[i]))
        # If flyped crossing begins at index 1 (including wraparound)
        elif (original_pd[i][1] < original_pd[i][3]) or (original_pd[i][1] == num_edges - 1):
            flyped_pd.append([(original_pd[i][1] - 1) % num_edges, (original_pd[i][0] - 1) % num_edges, (original_pd[i][3] - 1) % num_edges, (original_pd[i][2] - 1) % num_edges])
        # If flyped crossing begins at index 3
        else:
            flyped_pd.append([(original_pd[i][3] - 1) % num_edges, (original_pd[i][2] - 1) % num_edges, (original_pd[i][1] - 1) % num_edges, (original_pd[i][0] - 1) % num_edges])
    # Adding in flipped initial crossing; understrand depends on if in1/out1 are both top/bottom strands or if one is top, other is # bottom
    if out1 % 2 == 1:
        if original_pd[0][1] < original_pd[0][3] or original_pd[0][1] == num_edges - 1: 
            flyped_pd.append([(out1 - 1) % num_edges, (out2 - 1) % num_edges, out1, out2])
        else:
            flyped_pd.append([(out1 - 1) % num_edges, out2, out1, (out2 - 1) % num_edges])
    else:
        if original_pd[0][1] < original_pd[0][3] or original_pd[0][1] == num_edges - 1:
            flyped_pd.append([(out2 - 1) % num_edges, (out1 - 1) % num_edges, out2, out1])
        else:
            flyped_pd.append([(out2 - 1) % num_edges, out1, out2, (out1 - 1) % num_edges])
    flyped_pd.sort()
    return Link(flyped_pd)