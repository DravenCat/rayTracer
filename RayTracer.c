/*
  CSC D18 - RayTracer code.

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO" - remember to check what
  functionality is actually needed for the corresponding
  assignment!

  Last updated: Aug. 2017   - F.J.E.
*/

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name:	Dezhi Ren
* 2) Student Name:		
*
* 1) Student number: 1005736795
* 2) Student number:
* 
* 1) UtorID rendezhi
* 2) UtorID
* 
* We hereby certify that the work contained here is our own
*  REN,DEZHI
* ____________________             _____________________
* (sign with your name)            (sign with your name)
********************************************************************************/

#include "utils.h"    // <-- This includes RayTracer.h

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;

void buildScene(void) {
#include "buildscene.c"        // <-- Import the scene definition!
}

void calculatePhongModel(struct point3D ls_ray_d, struct ray3D *ray, struct pointLS *i, struct object3D *obj,
                         struct point3D *n, double R, double G, double B,
                         struct colourRGB *tmp_col) {
    // diffusion component rd*Id*max(0,n*s)
    normalize(&ls_ray_d);
    double diffusion_param =  max(0, fabs(dot(n, &ls_ray_d)));
    if (!obj->frontAndBack) diffusion_param =  max(0, dot(n, &ls_ray_d));
    tmp_col->R += R * obj->alb.rd * i->col.R * diffusion_param;
    tmp_col->G += G * obj->alb.rd * i->col.G * diffusion_param;
    tmp_col->B += B * obj->alb.rd * i->col.B * diffusion_param;

    //specular component rs*Is*max(0,c*m)^a
    // get m
    struct point3D ls_ray_d_inv;
    double dot_m = 2 * dot(n, &ls_ray_d);
    ls_ray_d_inv.px = -ls_ray_d.px + dot_m *n->px;
    ls_ray_d_inv.py = -ls_ray_d.py + dot_m *n->py;
    ls_ray_d_inv.pz = -ls_ray_d.pz + dot_m *n->pz;
    ls_ray_d_inv.pw = 1;
    normalize(&ls_ray_d_inv);

    // get c
    struct point3D ray_inv_d;
    ray_inv_d.px = - ray->d.px;
    ray_inv_d.py = - ray->d.py;
    ray_inv_d.pz = - ray->d.pz;
    ray_inv_d.pw = 1;
    normalize(&ray_inv_d);

    // get max(0, c*m)^a
    double specular_param = pow(fabs(dot(&ray_inv_d, &ls_ray_d_inv)), obj->shinyness);
    if (!obj->frontAndBack) specular_param = pow(max(0, dot(&ray_inv_d, &ls_ray_d_inv)), obj->shinyness);
    tmp_col->R += obj->alb.rs * i->col.R * specular_param;
    tmp_col->G += obj->alb.rs * i->col.G * specular_param;
    tmp_col->B += obj->alb.rs * i->col.B * specular_param;
}

void
rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b,
        struct colourRGB *col) {
    // This function implements the shading model as described in lecture. It takes
    // - A pointer to the first object intersected by the ray (to get the colour properties)
    // - The coordinates of the intersection point (in world coordinates)
    // - The normal at the point
    // - The ray (needed to determine the reflection direction to use for the global component, as well as for
    //   the Phong specular component)
    // - The current racursion depth
    // - The (a,b) texture coordinates (meaningless unless texture is enabled)
    //
    // Returns:
    // - The colour for this ray (using the col pointer)
    //

    struct colourRGB tmp_col;    // Accumulator for colour components
    double R, G, B;            // Colour for the object in R G and B

    // This will hold the colour as we process all the components of
    // the Phong illumination model
    tmp_col.R = 0;
    tmp_col.G = 0;
    tmp_col.B = 0;

    if (obj->texImg == NULL)        // Not textured, use object colour
    {
        R = obj->col.R;
        G = obj->col.G;
        B = obj->col.B;
    } else {
        // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
        // for the object. Note that we will use textures also for Photon Mapping.
        obj->textureMap(obj->texImg, a, b, &R, &G, &B);
    }

    //////////////////////////////////////////////////////////////
    // TO DO: Implement this function. Refer to the notes for
    // details about the shading model.
    //////////////////////////////////////////////////////////////

    // Be sure to update 'col' with the final colour computed here!
    double lambda_shade = -1.0;
    struct object3D *shade_obj;
    struct point3D shade_p, shade_n;
    double shade_a, shade_b;

    // Local component
    for (struct pointLS *i = light_list; i != NULL; i = i->next) {
        // initialize the ray from the light source
        struct ray3D light_src_ray;
        struct point3D d;
        memcpy(&d, &i->p0, sizeof(struct point3D));
        subVectors(p, &d);
        initRay(&light_src_ray, p, &d);

        // ambient component
        tmp_col.R += obj->alb.ra * R;
        tmp_col.G += obj->alb.ra * G;
        tmp_col.B += obj->alb.ra * B;

        // check whether in the shadow
        findFirstHit(&light_src_ray, &lambda_shade, obj, &shade_obj,
                     &shade_p, &shade_n, &shade_a, &shade_b);

        // lambda < 0 or > 1 means not in shadow
        if (lambda_shade + 1e-6 < 0 || lambda_shade > 1e-6 + 1) {
            // calculate phong model
            calculatePhongModel(light_src_ray.d, ray, i, obj, n, R, G, B, &tmp_col);
        }
    }

    // Global component
    if (depth < MAX_DEPTH) {

        // secondary reflection
        if (obj->alb.rg != 0){
            struct colourRGB mirror_result;
            // Get the ray in mirror direction
            struct ray3D mirror_ray;
            double m_dot = -2 * dot(&ray->d, n);
            mirror_ray.d.px = ray->d.px + m_dot * n->px;
            mirror_ray.d.py = ray->d.py + m_dot * n->py;
            mirror_ray.d.pz = ray->d.pz + m_dot * n->pz;
            mirror_ray.d.pw = 1;
            normalize(&mirror_ray.d);
            initRay(&mirror_ray, p, &mirror_ray.d);

            // Trace the ray and add its result
            rayTrace(&mirror_ray, depth + 1, &mirror_result, obj);
            tmp_col.R += mirror_result.R * obj->alb.rg;
            tmp_col.G += mirror_result.G * obj->alb.rg;
            tmp_col.B += mirror_result.B * obj->alb.rg;
        }
        // refraction
    }
    col->R = tmp_col.R;
    col->G = tmp_col.G;
    col->B = tmp_col.B;
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p,
                  struct point3D *n, double *a, double *b) {
    // Find the closest intersection between the ray and any objects in the scene.
    // Inputs:
    //   *ray    -  A pointer to the ray being traced
    //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
    //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
    //              projection
    // Outputs:
    //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
    //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
    //              this ray (this is required so you can do the shading)
    //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
    //   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
    //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

    /////////////////////////////////////////////////////////////
    // TO DO: Implement this function. See the notes for
    // reference of what to do in here
    /////////////////////////////////////////////////////////////
    double lambda_min = -1.0;
    *lambda = -1.0;

    // traverse all the object
    for (struct object3D *i = object_list; i != NULL; i = i->next) {
        if (i != Os){
            double obj_lambda = -1;

            //find intersection
            (i->intersect)(i, ray, &obj_lambda, p, n, a, b);

            // replace first hit with the smallest positive lambda
            if (((lambda_min + 1e-6 <= 0) || (obj_lambda + 1e-6 < lambda_min)) && obj_lambda > 1e-6) {
                lambda_min = obj_lambda;
                *(lambda) = lambda_min;
                *(obj) = i;
            }
        }
    }
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os) {
    // Trace one ray through the scene.
    //
    // Parameters:
    //   *ray   -  A pointer to the ray being traced
    //   depth  -  Current recursion depth for recursive raytracing
    //   *col   - Pointer to an RGB colour structure so you can return the object colour
    //            at the intersection point of this ray with the closest scene object.
    //   *Os    - 'Object source' is a pointer to the object from which the ray
    //            originates so you can discard self-intersections due to numerical
    //            errors. NULL for rays originating from the center of projection.

    double lambda;        // Lambda at intersection
    double a, b;        // Texture coordinates
    struct object3D *obj;    // Pointer to object at intersection
    struct point3D p;    // Intersection point
    struct point3D n;    // Normal at intersection
    struct colourRGB I;    // Colour returned by shading function

    if (depth > MAX_DEPTH)    // Max recursion depth reached. Return invalid colour.
    {
        col->R = -1;
        col->G = -1;
        col->B = -1;
        return;
    }

    ///////////////////////////////////////////////////////
    // TO DO: Complete this function. Refer to the notes
    // if you are unsure what to do here.
    ///////////////////////////////////////////////////////

    findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
    if (lambda > 1e-6) { // if hit an object
        rtShade(obj, &p, &n, ray, depth, a, b, &I);
        // colRGB in  [0, 1]
        col->R = max(0, I.R);
        col->G = max(0, I.G);
        col->B = max(0, I.B);
        col->R = min(1, I.R);
        col->G = min(1, I.G);
        col->B = min(1, I.B);
    } else { // hit nothing
        col->R = 0;
        col->G = 0;
        col->B = 0;
    }
}

int main(int argc, char *argv[]) {
    // Main function for the raytracer. Parses input parameters,
    // sets up the initial blank image, and calls the functions
    // that set up the scene and do the raytracing.
    struct image *im;    // Will hold the raytraced image
    struct view *cam;    // Camera and view for this scene
    int sx;        // Size of the raytraced image
    int antialiasing;    // Flag to determine whether antialiaing is enabled or disabled
    char output_name[1024];    // Name of the output file for the raytraced .ppm image
    struct point3D e;        // Camera view parameters 'e', 'g', and 'up'
    struct point3D g;
    struct point3D up;
    double du, dv;            // Increase along u and v directions for pixel coordinates
    struct point3D pc, d;        // Point structures to keep the coordinates of a pixel and
    // the direction or a ray
    struct ray3D ray;        // Structure to keep the ray from e to a pixel
    struct colourRGB col;        // Return colour for raytraced pixels
    struct colourRGB background;   // Background colour
    int i, j;            // Counters for pixel coordinates
    unsigned char *rgbIm;

    if (argc < 5) {
        fprintf(stderr, "RayTracer: Can not parse input parameters\n");
        fprintf(stderr, "USAGE: RayTracer size rec_depth antialias output_name\n");
        fprintf(stderr, "   size = Image size (both along x and y)\n");
        fprintf(stderr, "   rec_depth = Recursion depth\n");
        fprintf(stderr, "   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
        fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
        exit(0);
    }
    sx = atoi(argv[1]);
    MAX_DEPTH = atoi(argv[2]);
    if (atoi(argv[3]) == 0) antialiasing = 0; else antialiasing = 1;
    strcpy(&output_name[0], argv[4]);

    fprintf(stderr, "Rendering image at %d x %d\n", sx, sx);
    fprintf(stderr, "Recursion depth = %d\n", MAX_DEPTH);
    if (!antialiasing) fprintf(stderr, "Antialising is off\n");
    else fprintf(stderr, "Antialising is on\n");
    fprintf(stderr, "Output file name: %s\n", output_name);

    object_list = NULL;
    light_list = NULL;
    texture_list = NULL;

    // Allocate memory for the new image
    im = newImage(sx, sx);
    if (!im) {
        fprintf(stderr, "Unable to allocate memory for raytraced image\n");
        exit(0);
    } else rgbIm = (unsigned char *) im->rgbdata;

    ///////////////////////////////////////////////////
    // TO DO: You will need to implement several of the
    //        functions below. For Assignment 2, you can use
    //        the simple scene already provided. But
    //        for Assignment 3 you need to create your own
    //        *interesting* scene.
    ///////////////////////////////////////////////////
    buildScene();        // Create a scene. This defines all the
    // objects in the world of the raytracer

    //////////////////////////////////////////
    // TO DO: For Assignment 2 you can use the setup
    //        already provided here. For Assignment 3
    //        you may want to move the camera
    //        and change the view parameters
    //        to suit your scene.
    //////////////////////////////////////////

    // Mind the homogeneous coordinate w of all vectors below. DO NOT
    // forget to set it to 1, or you'll get junk out of the
    // geometric transformations later on.

    // Camera center is at (0,0,-1)
    e.px = 0;
    e.py = 0;
    e.pz = -1;
    e.pw = 1;

    // To define the gaze vector, we choose a point 'pc' in the scene that
    // the camera is looking at, and do the vector subtraction pc-e.
    // Here we set up the camera to be looking at the origin.
    g.px = 0 - e.px;
    g.py = 0 - e.py;
    g.pz = 0 - e.pz;
    g.pw = 1;
    // In this case, the camera is looking along the world Z axis, so
    // vector w should end up being [0, 0, -1]

    // Define the 'up' vector to be the Y axis
    up.px = 0;
    up.py = 1;
    up.pz = 0;
    up.pw = 1;

    // Set up view with given the above vectors, a 4x4 window,
    // and a focal length of -1 (why? where is the image plane?)
    // Note that the top-left corner of the window is at (-2, 2)
    // in camera coordinates.
    cam = setupView(&e, &g, &up, -1, -2, 2, 4);

    if (cam == NULL) {
        fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
        cleanup(object_list, light_list, texture_list);
        deleteImage(im);
        exit(0);
    }

    // Set up background colour here
    background.R = 0;
    background.G = 0;
    background.B = 0;

    // Do the raytracing
    //////////////////////////////////////////////////////
    // TO DO: You will need code here to do the raytracing
    //        for each pixel in the image. Refer to the
    //        lecture notes, in particular, to the
    //        raytracing pseudocode, for details on what
    //        to do here. Make sure you undersand the
    //        overall procedure of raytracing for a single
    //        pixel.
    //////////////////////////////////////////////////////
    du = cam->wsize / (sx - 1);        // du and dv. In the notes in terms of wl and wr, wt and wb,
    dv = -cam->wsize / (sx - 1);        // here we use wl, wt, and wsize. du=dv since the image is
    // and dv is negative since y increases downward in pixel
    // coordinates and upward in camera coordinates.

    fprintf(stderr, "View parameters:\n");
    fprintf(stderr, "Left=%f, Top=%f, Width=%f, f=%f\n", cam->wl, cam->wt, cam->wsize, cam->f);
    fprintf(stderr, "Camera to world conversion matrix (make sure it makes sense!):\n");
    printmatrix(cam->C2W);
    fprintf(stderr, "World to camera conversion matrix:\n");
    printmatrix(cam->W2C);
    fprintf(stderr, "\n");

    fprintf(stderr, "Rendering row: ");
    for (j = 0; j < sx; j++)        // For each of the pixels in the image
    {
        fprintf(stderr, "%d/%d, ", j, sx);
        for (i = 0; i < sx; i++) {
            ///////////////////////////////////////////////////////////////////
            // TO DO - complete the code that should be in this loop to do the
            //         raytracing!
            ///////////////////////////////////////////////////////////////////
            // Calculate pc in camera coordinate
            pc.px = cam->wl + i * du;
            pc.py = cam->wt + j * dv;
            pc.pz = cam->f;
            pc.pw = 1;

            // Convert pc into pw and calculate the direction of ray
            matVecMult(cam->C2W, &pc); // pw = M_cw' * pc
            memcpy(&d, &pc, sizeof(struct point3D));
            subVectors(&e, &d); // d = pw - e
            normalize(&d);

            // Initialize the ray
            initRay(&ray, &cam->e, &d);
            normalize(&ray.d);

            // Trace the ray of the pixel
            memcpy(&col, &background, sizeof(struct colourRGB));
            rayTrace(&ray, 1, &col, NULL);

            // Set the final colorRGB of the pixel
            *(rgbIm + 3 * (sx * j + i)) = (unsigned char) (col.R * 255.0); // Red
            *(rgbIm + 3 * (sx * j + i) + 1) = (unsigned char) (col.G * 255.0); // Green
            *(rgbIm + 3 * (sx * j + i) + 2) = (unsigned char) (col.B * 255.0); // Blue

        } // end for i
    } // end for j

    fprintf(stderr, "\nDone!\n");

    // Output rendered image
    imageOutput(im, output_name);

    // Exit section. Clean up and return.
    cleanup(object_list, light_list, texture_list);        // Object, light, and texture lists
    deleteImage(im);                    // Rendered image
    free(cam);                        // camera view
    exit(0);
}

