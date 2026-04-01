#include <gtest/gtest.h>
#include <DSTmath/DSTgeo.h>

using namespace DST::Math;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Build a unit square with corners (0,0)-(1,0)-(1,1)-(0,1). */
static polygon2D makeUnitSquare()
{
    polygon2D p;
    p.AddPoint(point({0., 0.}));
    p.AddPoint(point({1., 0.}));
    p.AddPoint(point({1., 1.}));
    p.AddPoint(point({0., 1.}));
    return p;
}

/** Build a 2×2 square with corners (0,0)-(2,0)-(2,2)-(0,2). */
static polygon2D make2x2Square()
{
    polygon2D p;
    p.AddPoint(point({0., 0.}));
    p.AddPoint(point({2., 0.}));
    p.AddPoint(point({2., 2.}));
    p.AddPoint(point({0., 2.}));
    return p;
}

// ---------------------------------------------------------------------------
#pragma region - ctor tests

TEST(polygon2D, default_ctor)
{
    polygon2D p;
    EXPECT_EQ(p.size(), 0u);
}

TEST(polygon2D, copy_ctor)
{
    polygon2D sq = makeUnitSquare();
    EXPECT_EQ(sq.size(), 4u);

    polygon2D cp(sq);
    EXPECT_EQ(cp.size(), sq.size());

    for (size_t i = 0; i < sq.size(); i++)
    {
        EXPECT_NEAR(cp[i].X(), sq[i].X(), point::precision);
        EXPECT_NEAR(cp[i].Y(), sq[i].Y(), point::precision);
    }
}

TEST(polygon2D, copy_ctor_empty)
{
    polygon2D empty;
    polygon2D cp(empty);
    EXPECT_EQ(cp.size(), 0u);
}

#pragma endregion

// ---------------------------------------------------------------------------
#pragma region - AddPoint tests

TEST(polygon2D, AddPoint_basic)
{
    polygon2D p;
    EXPECT_EQ(p.size(), 0u);

    p.AddPoint(point({1., 2.}));
    EXPECT_EQ(p.size(), 1u);
    EXPECT_NEAR(p[0].X(), 1., point::precision);
    EXPECT_NEAR(p[0].Y(), 2., point::precision);

    p.AddPoint(point({3., 4.}));
    EXPECT_EQ(p.size(), 2u);

    p.AddPoint(point({5., 6.}));
    EXPECT_EQ(p.size(), 3u);
}

TEST(polygon2D, AddPoint_rejects_duplicate)
{
    polygon2D p;
    p.AddPoint(point({1., 2.}));
    p.AddPoint(point({3., 4.}));
    EXPECT_EQ(p.size(), 2u);

    // Adding exact duplicates must not grow the polygon.
    p.AddPoint(point({1., 2.}));
    EXPECT_EQ(p.size(), 2u);

    p.AddPoint(point({3., 4.}));
    EXPECT_EQ(p.size(), 2u);
}

TEST(polygon2D, AddPoint_stores_only_XY)
{
    // A 3D point: only X and Y should be stored.
    polygon2D p;
    p.AddPoint(point({7., 8., 99.}));
    EXPECT_EQ(p.size(), 1u);
    EXPECT_NEAR(p[0].X(), 7., point::precision);
    EXPECT_NEAR(p[0].Y(), 8., point::precision);
    // The stored vertex is 2D, so it has no Z component.
    EXPECT_EQ(p[0].size(), 2u);
}

#pragma endregion

// ---------------------------------------------------------------------------
#pragma region - Area tests

TEST(polygon2D, Area_empty)
{
    polygon2D p;
    EXPECT_NEAR(p.Area(), 0., point::precision);
}

TEST(polygon2D, Area_too_few_vertices)
{
    polygon2D p;
    p.AddPoint(point({0., 0.}));
    p.AddPoint(point({1., 0.}));
    EXPECT_NEAR(p.Area(), 0., point::precision);
}

TEST(polygon2D, Area_unit_square)
{
    polygon2D sq = makeUnitSquare();
    EXPECT_NEAR(sq.Area(), 1.0, 1e-10);
}

TEST(polygon2D, Area_2x2_square)
{
    polygon2D sq = make2x2Square();
    EXPECT_NEAR(sq.Area(), 4.0, 1e-10);
}

TEST(polygon2D, Area_rectangle)
{
    // 2 × 3 rectangle → area = 6
    polygon2D rect;
    rect.AddPoint(point({0., 0.}));
    rect.AddPoint(point({2., 0.}));
    rect.AddPoint(point({2., 3.}));
    rect.AddPoint(point({0., 3.}));
    EXPECT_NEAR(rect.Area(), 6.0, 1e-10);
}

TEST(polygon2D, Area_right_triangle)
{
    // Right triangle (0,0)-(4,0)-(0,3) → area = 0.5×4×3 = 6
    polygon2D tri;
    tri.AddPoint(point({0., 0.}));
    tri.AddPoint(point({4., 0.}));
    tri.AddPoint(point({0., 3.}));
    EXPECT_NEAR(tri.Area(), 6.0, 1e-10);
}

TEST(polygon2D, Area_inverted_vertex_order)
{
    // Adding vertices counter-clockwise should give the same area after reordering.
    polygon2D sq;
    sq.AddPoint(point({0., 1.}));
    sq.AddPoint(point({1., 1.}));
    sq.AddPoint(point({1., 0.}));
    sq.AddPoint(point({0., 0.}));
    EXPECT_NEAR(sq.Area(), 1.0, 1e-10);
}

#pragma endregion

// ---------------------------------------------------------------------------
#pragma region - IsPointInside tests

TEST(polygon2D, IsPointInside_empty_polygon)
{
    polygon2D p;
    EXPECT_FALSE(p.IsPointInside(point({0., 0.})));
}

TEST(polygon2D, IsPointInside_two_vertices)
{
    polygon2D p;
    p.AddPoint(point({0., 0.}));
    p.AddPoint(point({1., 0.}));
    EXPECT_FALSE(p.IsPointInside(point({0.5, 0.})));
}

TEST(polygon2D, IsPointInside_center)
{
    polygon2D sq = makeUnitSquare();
    EXPECT_TRUE(sq.IsPointInside(point({0.5, 0.5})));
}

TEST(polygon2D, IsPointInside_clearly_outside)
{
    polygon2D sq = makeUnitSquare();
    EXPECT_FALSE(sq.IsPointInside(point({ 2.0,  0.5})));
    EXPECT_FALSE(sq.IsPointInside(point({-1.0,  0.5})));
    EXPECT_FALSE(sq.IsPointInside(point({ 0.5, -1.0})));
    EXPECT_FALSE(sq.IsPointInside(point({ 0.5,  2.0})));
}

TEST(polygon2D, IsPointInside_on_edges)
{
    polygon2D sq = makeUnitSquare();
    EXPECT_TRUE(sq.IsPointInside(point({0.5, 0.0})));  // bottom edge
    EXPECT_TRUE(sq.IsPointInside(point({0.5, 1.0})));  // top edge
    EXPECT_TRUE(sq.IsPointInside(point({0.0, 0.5})));  // left edge
    EXPECT_TRUE(sq.IsPointInside(point({1.0, 0.5})));  // right edge
}

TEST(polygon2D, IsPointInside_on_corners)
{
    polygon2D sq = makeUnitSquare();
    EXPECT_TRUE(sq.IsPointInside(point({0., 0.})));
    EXPECT_TRUE(sq.IsPointInside(point({1., 0.})));
    EXPECT_TRUE(sq.IsPointInside(point({1., 1.})));
    EXPECT_TRUE(sq.IsPointInside(point({0., 1.})));
}

TEST(polygon2D, IsPointInside_triangle)
{
    polygon2D tri;
    tri.AddPoint(point({0., 0.}));
    tri.AddPoint(point({4., 0.}));
    tri.AddPoint(point({0., 4.}));

    EXPECT_TRUE (tri.IsPointInside(point({1.,  1.})));  // inside
    EXPECT_FALSE(tri.IsPointInside(point({3.,  3.})));  // outside (beyond hypotenuse)
    EXPECT_TRUE (tri.IsPointInside(point({2.,  0.})));  // on base edge
    EXPECT_TRUE (tri.IsPointInside(point({0.,  2.})));  // on left edge
}

#pragma endregion

// ---------------------------------------------------------------------------
#pragma region - GetIntersectionWithPolygons tests

TEST(polygon2D, Intersection_throws_on_small_polygon)
{
    polygon2D small;
    small.AddPoint(point({0., 0.}));
    small.AddPoint(point({1., 0.}));

    polygon2D valid = makeUnitSquare();

    EXPECT_ANY_THROW(small.GetIntersectionWithPolygons(valid));
    EXPECT_ANY_THROW(valid.GetIntersectionWithPolygons(small));
}

TEST(polygon2D, Intersection_non_overlapping)
{
    // Two squares far apart — intersection must be empty (area == 0).
    polygon2D sq1 = makeUnitSquare();

    polygon2D sq2;
    sq2.AddPoint(point({5., 5.}));
    sq2.AddPoint(point({6., 5.}));
    sq2.AddPoint(point({6., 6.}));
    sq2.AddPoint(point({5., 6.}));

    polygon2D inter = sq1.GetIntersectionWithPolygons(sq2);
    EXPECT_NEAR(inter.Area(), 0., 1e-10);
}

TEST(polygon2D, Intersection_identical_squares)
{
    // Intersection of a square with itself must equal the square area.
    polygon2D sq1 = makeUnitSquare();
    polygon2D sq2 = makeUnitSquare();

    polygon2D inter = sq1.GetIntersectionWithPolygons(sq2);
    EXPECT_NEAR(inter.Area(), 1.0, 1e-10);
}

TEST(polygon2D, Intersection_half_overlap)
{
    // sq1: (0,0)-(2,0)-(2,2)-(0,2), area=4
    // sq2: (1,0)-(3,0)-(3,2)-(1,2), area=4
    // Overlap: (1,0)-(2,0)-(2,2)-(1,2), area=2
    polygon2D sq1 = make2x2Square();

    polygon2D sq2;
    sq2.AddPoint(point({1., 0.}));
    sq2.AddPoint(point({3., 0.}));
    sq2.AddPoint(point({3., 2.}));
    sq2.AddPoint(point({1., 2.}));

    polygon2D inter = sq1.GetIntersectionWithPolygons(sq2);
    EXPECT_NEAR(inter.Area(), 2.0, 1e-10);
}

TEST(polygon2D, Intersection_quarter_overlap)
{
    // sq1: (0,0)-(2,0)-(2,2)-(0,2), area=4
    // sq2: (1,1)-(3,1)-(3,3)-(1,3), area=4
    // Overlap: (1,1)-(2,1)-(2,2)-(1,2), area=1
    polygon2D sq1 = make2x2Square();

    polygon2D sq2;
    sq2.AddPoint(point({1., 1.}));
    sq2.AddPoint(point({3., 1.}));
    sq2.AddPoint(point({3., 3.}));
    sq2.AddPoint(point({1., 3.}));

    polygon2D inter = sq1.GetIntersectionWithPolygons(sq2);
    EXPECT_EQ(inter.size(), 4u);
    EXPECT_NEAR(inter.Area(), 1.0, 1e-10);
}

TEST(polygon2D, Intersection_one_inside_other)
{
    // Small square fully inside big square — intersection equals small square.
    polygon2D big = make2x2Square();

    polygon2D small = makeUnitSquare();   // (0,0)-(1,0)-(1,1)-(0,1)

    polygon2D inter = big.GetIntersectionWithPolygons(small);
    EXPECT_NEAR(inter.Area(), 1.0, 1e-10);
}

#pragma endregion

// ---------------------------------------------------------------------------
#pragma region - Dump tests

TEST(polygon2D, Dump_empty)
{
    polygon2D p;
    std::string s = p.Dump("empty_poly");
    EXPECT_FALSE(s.empty());
}

TEST(polygon2D, Dump_non_empty)
{
    polygon2D sq = makeUnitSquare();
    std::string s = sq.Dump("square");
    EXPECT_FALSE(s.empty());
}

TEST(polygon2D, DumpCorner_empty)
{
    polygon2D p;
    std::string s = p.DumpCorner();
    EXPECT_FALSE(s.empty());
}

TEST(polygon2D, DumpCorner_non_empty)
{
    polygon2D sq = makeUnitSquare();
    std::string s = sq.DumpCorner();
    EXPECT_FALSE(s.empty());
    // The first vertex must appear twice (open and closed path).
    size_t first = s.find('(');
    size_t last  = s.rfind('(');
    EXPECT_NE(first, last);
}

#pragma endregion
