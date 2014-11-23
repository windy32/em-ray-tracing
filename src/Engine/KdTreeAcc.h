#ifndef KD_TREE_ACC_H
#define KD_TREE_ACC_H

#include "Accelerator.h"

class KdTreeAcc : public Accelerator
{
private:
    enum Axes { XAxis, YAxis, ZAxis, NoAxis }; // "NoAxis" denotes a leaf
    enum KdEventType { End, Planar, Start };

    struct KdNode
    {
        std::vector<Geometry *> list; // list of enclosed objects
        KdNode *left;  // pointer to the left child
        KdNode *right; // pointer to the right child
        Axes axis;         // orientation of the splitting plane
        double splitPlane;  // position of the splitting plane

        Point min;
        Point max;
    };
    KdNode *root;

    struct StackElem
    {
        KdNode *node;  // pointer of far child
        double t;           // the entry / exit signed distance
        Point pb;          // the coordinates of entry / exit point
        int prev;          // the pointer to the previous stack item
    };

public: // should be exposed to the compare functions for sorting
    struct KdEvent
    {
        Geometry *triangle;
        double position;
        KdEventType type;

        KdEvent(Geometry *triangle, double position, KdEventType type) : 
            triangle(triangle), position(position), type(type) {}
    };

private:
    void buildKdTree(KdNode *node, std::vector<Geometry *> &list, int depth, int &numLeaves, int &leafElements);
    void deleteTree(KdNode *node);
    double split(KdNode *node, int axis, std::vector<Geometry *> &list);
    double splitSAH(KdNode *node, std::vector<Geometry *> &list, int &bestAxis, double &minSAH);

public:
    KdTreeAcc(std::vector<Geometry *> *scene) : Accelerator(scene) {}
    ~KdTreeAcc();
    virtual void init();
    virtual IntersectResult intersect(Ray &ray, std::vector<RxIntersection> &rxPoints);
};

#endif
