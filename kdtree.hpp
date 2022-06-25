/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <iostream>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim> &first,
                                const Point<Dim> &second, int curDim) const
{
  if (first[curDim] == second[curDim])
  {
    return first<second;
  }

  return first[curDim] < second[curDim];
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim> &target,
                                const Point<Dim> &currentBest,
                                const Point<Dim> &potential) const
{
  /**
   * @todo Implement this function!
   */
  double curDist = 0.0;
  double potDist = 0.0;
  for (int i = 0; i < Dim; i++)
  {
    curDist += std::pow(currentBest[i] - target[i], 2.0);
    potDist += std::pow(potential[i] - target[i], 2.0);
  }
  if (curDist==potDist)
  {
    // std::cout<<"tie: "<<curDist<<" "<<potDist<<std::endl;
    // std::cout<<"tie: "<<currentBest<<" "<<potential<<" "<<target<<std::endl;
    // if(potential<currentBest) {
    // std::cout<< potential <<std::endl;
    // } else {
    //   std::cout<< currentBest <<std::endl;
    // }
    return potential<currentBest;
  }
  return potDist < curDist;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>> &newPoints)
{
  /**
   * @todo Implement this function!
   * root = median node
   * partition vector
   * root->left = median node of left half
   * root->right = median node of right half
   */
  vector<Point<Dim>> pts = newPoints;
  root = make_tree(pts, 0, pts.size() - 1, 0, root);
}
template <int Dim>
typename KDTree<Dim>::KDTreeNode *KDTree<Dim>::make_tree(vector<Point<Dim>> &pts, int start, int end, int d, KDTreeNode *root)
{
  if (start <= end)
  {
    int mid = (end + start) / 2;
    Point<Dim> midPoint = select(pts, start, end, d, mid);
    root = new KDTreeNode(midPoint);

    root->left = make_tree(pts, start, mid - 1, (d + 1) % Dim, root->left);
    root->right = make_tree(pts, mid + 1, end, (d + 1) % Dim, root->right);
    return root;
  }
  else
  {
    return nullptr;
  }
}

template <int Dim>
Point<Dim> KDTree<Dim>::select(vector<Point<Dim>> &pts, int start, int end, int d, int k)
{
  int pivotIndex = 0;
  while (true)
  {
    if (start == end)
    { // If the list contains only one element,
      return pts[start];
    }
    int pivotIndex = k;
    pivotIndex = partition(pts, start, end, d, pivotIndex);

    if (pivotIndex == k)
    {
      return pts[pivotIndex];
    }
    else if (k < pivotIndex)
    {
      end = pivotIndex - 1;
    }
    else
    {
      start = pivotIndex + 1;
    }
  }
  return pts[k];
}

template <int Dim>
int KDTree<Dim>::partition(vector<Point<Dim>> &pts, int start, int end, int d, int k)
{
  Point<Dim> pivotValue = pts[k];
  iter_swap(pts.begin() + k, pts.begin() + end);
  int temp = start;
  for (int i = start; i < end; i++)
  {
    if (smallerDimVal(pts[i], pivotValue, d))
    {
      iter_swap(pts.begin() + temp, pts.begin() + i);
      temp++;
    }
  }
  iter_swap(pts.begin() + temp, pts.begin() + end);
  return temp;
}

template <int Dim>
void KDTree<Dim>::_delete()
{
  if (root)
  {
    delete root;
  }
}

template <int Dim>
void KDTree<Dim>::_copy(const KDTree<Dim> &other)
{
  root = new KDTreeNode(other.root);
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim> &other)
{
  /**
   * @todo Implement this function!
   */
  _copy(other);
}

template <int Dim>
const KDTree<Dim> &KDTree<Dim>::operator=(const KDTree<Dim> &rhs)
{

  /**
   * @todo Implement this function!
   */
  if (this != &rhs)
  {
    _delete();
    _copy(rhs);
  }

  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree()
{
  /**
   * @todo Implement this function!
   */
  _delete();
  root = nullptr;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(KDTreeNode *curRoot, Point<Dim> query, int d) const
{
  Point<Dim> nearest;
  // printTree();
  // std::cout<<curRoot->point<<std::endl;
    if(!curRoot->left && !curRoot->right ) return curRoot->point;

    if (smallerDimVal(query, curRoot->point, d))
    {
      if (curRoot->left)
      {
        nearest = findNearestNeighbor(curRoot->left, query, (d + 1) % Dim);
      }
      else
      {
        nearest = findNearestNeighbor(curRoot->right, query, (d + 1) % Dim);
      }
    }
    else
    {
      if (curRoot->right)
      {
        nearest = findNearestNeighbor(curRoot->right, query, (d + 1) % Dim);
      }
      else
      { 
        nearest = findNearestNeighbor(curRoot->left, query, (d + 1) % Dim);
      }
    }

    if (shouldReplace(query, nearest, curRoot->point))
    {
      nearest = curRoot->point;
    }

    double radius = 0;
    for (int i = 0; i < Dim; i++)
    {
      radius += std::pow(nearest[i] - query[i], 2.0);
    }
    double splitDist = std::pow(curRoot->point[d] - query[d], 2.0);

    if (radius >= splitDist)
    {
      Point<Dim> tempNearest;
      if (smallerDimVal(query, curRoot->point, d) && curRoot->right)
      {
        tempNearest = findNearestNeighbor(curRoot->right, query, (d + 1) % Dim);
        if (shouldReplace(query, nearest, tempNearest))
      {
        nearest = tempNearest;
      else if (curRoot->left)
      {
        tempNearest = findNearestNeighbor(curRoot->left, query, (d + 1) % Dim);
        if (shouldReplace(query, nearest, tempNearest))
      {
        nearest = tempNearest;
      }
      }

      
    }
  
  return nearest;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim> &query) const
{
  /**
   * @todo Implement this function!
   */
  return findNearestNeighbor(root, query, 0);
}