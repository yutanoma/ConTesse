#include "tessellate_face.h"

namespace utils {

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXi,
           Eigen::MatrixXi, Eigen::MatrixXi, std::map<Point_2, int>,
           std::map<Segment_2, int, Segment2Comparator>>
tessellate_face(const Arrangement_with_history_2 &/*arr*/,
                const Face_const_handle &/*face*/)
{
  Eigen::MatrixXd V(0, 2);
  Eigen::MatrixXi F(0, 3);
  Eigen::MatrixXi EV(0, 2);
  Eigen::MatrixXi FE(0, 2);
  Eigen::MatrixXi EF(0, 3);
  std::map<Point_2, int> point_to_vid;
  std::map<Segment_2, int, Segment2Comparator> segment_to_eid;
  return {V, F, EV, FE, EF, point_to_vid, segment_to_eid};
}
}
