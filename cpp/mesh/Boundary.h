#include <vector>
#include <algorithm>

namespace mesh {
  class Boundary {
    public:
      int nfacets;
      std::vector<int> facets,//sorted boundary facets
        parentEls;// length of total # of facets. -1 if internal facet
                  // ielem of parent el if [iboun]
      std::vector<int> difference(const Boundary &otherBoundary) {
        /*
         * Return facets that are ours but not theirs.
         */
        std::vector<int> diff;
        std::set_difference(facets.begin(), facets.end(), otherBoundary.facets.begin(), otherBoundary.facets.end(),
                std::inserter(diff, diff.begin()));
        return diff;
      }
  };
}
