/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include "maptiles.h"
#include <iostream>
#include <map>
//#include "cs225/RGB_HSL.h"

using namespace std;

Point<3> convertToXYZ(LUVAPixel pixel) {
  return Point<3>(pixel.l, pixel.u, pixel.v);
}

MosaicCanvas *mapTiles(SourceImage const &theSource,
                       vector<TileImage> &theTiles) {
  /**
   * @todo Implement this function!
   */
  MosaicCanvas *out =
      new MosaicCanvas(theSource.getRows(), theSource.getColumns());

  map<Point<3>, TileImage*> map;

  // convert tile color to points vector and create map(avgColorPoint : Tile)
  vector<Point<3>> tileColors;

  for (unsigned int i = 0; i < theTiles.size(); i++) {
    Point<3> avgTileColor = Point<3>(theTiles[i].getAverageColor().l,
                                     theTiles[i].getAverageColor().u,
                                     theTiles[i].getAverageColor().v);
    map[avgTileColor] = &theTiles[i];
    tileColors.push_back(avgTileColor);
  }

  // create KDTree from Points vector
  KDTree<3> tileColorTree(tileColors);

  // use findNearestNeighbhor to find best tile to match region average
  for (int r = 0; r < theSource.getRows(); r++) {
    for (int c = 0; c < theSource.getColumns(); c++) {
      Point<3> avgRegionColor = Point<3>(theSource.getRegionColor(r, c).l,
                                         theSource.getRegionColor(r, c).u,
                                         theSource.getRegionColor(r, c).v);
      out->setTile(
          r, c, (map.at(tileColorTree.findNearestNeighbor(avgRegionColor))));
    }
  }

  return out;
}