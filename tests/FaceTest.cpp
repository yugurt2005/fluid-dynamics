#include <catch2/catch_test_macros.hpp>

#include "../src/Face.h"

TEST_CASE("Face: Normal") {
  SECTION("Horizontal") {
    Face face(-1, -1, {0, 1}, {0, 5});
    CHECK(face.normal == Vector2d(1, 0));
  }

  SECTION("Vertical") {
    Face face(-1, -1, {-3, 2}, {2, 2});
    CHECK(face.normal == Vector2d(0, -1));
  }

  SECTION("Diagonal") {
    Face face(-1, -1, {-1, -1}, {1, 1});
    CHECK(face.normal == Vector2d(1, -1).normalized());
  }
}

TEST_CASE("Face: Center") {
  SECTION("Horizontal") {
    Face face(-1, -1, {0, 1}, {0, 5});
    CHECK(face.center == Vector2d(0, 3));
  }

  SECTION("Vertical") {
    Face face(-1, -1, {-3, 2}, {2, 2});
    CHECK(face.center == Vector2d(-0.5, 2));
  }

  SECTION("Diagonal") {
    Face face(-1, -1, {-1, -1}, {1, 1});
    CHECK(face.center == Vector2d(0, 0));
  }
}

TEST_CASE("Face: Area") {
  SECTION("Horizontal") {
    Face face(-1, -1, {0, 1}, {0, 5});
    CHECK(face.area == 4);
  }

  SECTION("Vertical") {
    Face face(-1, -1, {-3, 2}, {2, 2});
    CHECK(face.area == 5);
  }

  SECTION("Diagonal") {
    Face face(-1, -1, {-1, -1}, {1, 1});
    CHECK(face.area == sqrt(8));
  }
}