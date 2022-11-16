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

void forwardPassTrace(struct ray3D *ray, int depth, struct object3D *Os, double R, double G, double B, int imgsize){
    if (depth > MAX_DEPTH) return;
    else {
        double tmp_lambda, tmp_a, tmp_b;
        struct object3D *first_hit;
        struct point3D tmp_p, tmp_n;
        findFirstHit(ray, &tmp_lambda, Os, &first_hit, &tmp_p, &tmp_n, &tmp_a, &tmp_b);

        if (tmp_lambda + 1e-6 < 0) { // hit nothing, stop tracing
            return;
        } else {
            // the ray hits a reflecting object.
            if (first_hit->alpha + 1e-6 >= 1 && first_hit->alb.rg > 1e-6) {
                double dot_product = dot(&ray->d, &tmp_n);

                struct point3D mirror_d;
                mirror_d.px = tmp_n.px * (-2) * dot_product + ray->d.px;
                mirror_d.py = tmp_n.py * (-2) * dot_product + ray->d.py;
                mirror_d.pz = tmp_n.pz * (-2) * dot_product + ray->d.pz;
                mirror_d.pw = 1;
                normalize(&mirror_d);

                ray = newRay(&tmp_p, &mirror_d);
                forwardPassTrace(ray, depth + 1, first_hit, first_hit->alb.rg * R, first_hit->alb.rg * G, first_hit->alb.rg * B, imgsize);
                free(ray);

            } else if (first_hit->alpha + 1e-6 < 1) {
                // the ray hits a refracting object.
                double r_idx1 = ray->insideOut ? 1.0 : first_hit->r_index;
                double r_idx2 = ray->insideOut ? first_hit->r_index : 1.0;
                struct point3D *op_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);
                double cos_theta1 = dot(op_d, &tmp_n);
                double sin_theta1 = sqrt(1 - pow(cos_theta1, 2));
                double sin_theta2 = (double) (r_idx1 / r_idx2) * sin_theta1;

                if (sin_theta2 < 1 && sin_theta2 > 0) {
                    double n21 = r_idx1 / r_idx2;
                    double dot_product = -dot(&tmp_n, &ray->d);
                    double tmp = sqrt(1 - n21 * n21 * (1 - dot_product * dot_product));

                    struct point3D *refract_d = newPoint(n21 * (dot_product * tmp_n.px + ray->d.px) - tmp * tmp_n.px,
                                                         n21 * (dot_product * tmp_n.py + ray->d.py) - tmp * tmp_n.py,
                                                         n21 * (dot_product * tmp_n.pz + ray->d.pz) - tmp * tmp_n.pz);
                    normalize(refract_d);
                    struct ray3D *refract_ray = newRay(&tmp_p, refract_d);
                    refract_ray->insideOut = 1 - ray->insideOut;

                    forwardPassTrace(refract_ray, depth + 1, first_hit,
                                     (1 - first_hit->alpha) * first_hit->col.R * R,
                                     (1 - first_hit->alpha) * first_hit->col.G * G,
                                     (1 - first_hit->alpha) * first_hit->col.B * B,
                                     imgsize);

                    free(refract_d);
                    free(refract_ray);
                }

                free(op_d);
            }

            if (first_hit->alb.rd && depth > 0) {
                // ray hits a diffuse surface & been reflected/refracted
                // store one photon at the location of the intersection point and stop tracing
                double *photon_rgb = (double *) first_hit->photonMap->rgbdata;

                int i = (int)tmp_a * imgsize;
                int j = (int)tmp_b * imgsize;

                *(photon_rgb + 3 * (i + imgsize * j)) += R;
                *(photon_rgb + 3 * (i + imgsize * j) + 1) += G;
                *(photon_rgb + 3 * (i + imgsize * j) + 2) += B;

//                photon_k =photon_k+1;
            }

        }
    }
}

void calculatePhongModel(struct point3D s, struct ray3D *ray, struct object3D *i, struct object3D *obj,
                         struct point3D *n, double R, double G, double B,
                         struct colourRGB *tmp_col) {
    // Diffuse term
    normalize(&s);
    double diffusion_param = fabs(dot(n, &s));
    if (!obj->frontAndBack && obj->alpha >= 1) diffusion_param = max(0, dot(n, &s));
    tmp_col->R += obj->alb.rd * R * i->col.R * diffusion_param;
    tmp_col->G += obj->alb.rd * G * i->col.G * diffusion_param;
    tmp_col->B += obj->alb.rd * B * i->col.B * diffusion_param;

    // Specular hightlights term

    struct point3D m;
    double dot_m = 2 * dot(n, &s);
    m.px = -s.px + dot_m * n->px;
    m.py = -s.py + dot_m * n->py;
    m.pz = -s.pz + dot_m * n->pz;
    m.pw = 1;
    normalize(&m);

    // get c
    struct point3D c;
    c.px = -ray->d.px;
    c.py = -ray->d.py;
    c.pz = -ray->d.pz;
    c.pw = 1;
    normalize(&c);
    double specular_param = pow(fabs(dot(&c, &m)), obj->shinyness);
    if (!obj->frontAndBack) specular_param = pow(max(0, dot(&c, &m)), obj->shinyness);
    tmp_col->R += obj->alb.rs * R * i->col.R * specular_param;
    tmp_col->G += obj->alb.rs * G * i->col.G * specular_param;
    tmp_col->B += obj->alb.rs * B * i->col.B * specular_param;
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
    double hasLS = 0;
    double lambda_shade = -1.0;
    struct object3D *shade_obj;
    struct point3D shade_p, shade_n;
    double shade_a, shade_b;


    struct point3D d;
    d.pw = -1;
    for (struct object3D *i = object_list; i != NULL; i = i->next) {
        if (i->isLightSource) {
            hasLS=1;
            int k = 0; // the number of unblocked rays
            int K = 5; // the number of samples
            for (int j = 0; j < K; j++) {
                // get a sample point on the object light source
                struct point3D p_obj;
                i->randomPoint(i, &p_obj.px, &p_obj.py, &p_obj.pz);

                // initialize the ray from the light source
                struct ray3D light_src_ray;
                memcpy(&d, &p_obj, sizeof(struct point3D));
                subVectors(p, &d);
                initRay(&light_src_ray, p, &d, 1);

                // check whether in the shadow
                findFirstHit(&light_src_ray, &lambda_shade, obj, &shade_obj,
                             &shade_p, &shade_n, &shade_a, &shade_b);

                // lambda < 0 or > 1 means not in shadow
                if (lambda_shade + 1e-6 < 0 || lambda_shade > 1e-6 + 1) {
                    k++;
                }
            }

            // exist some unblocked rays
            if (k > 0) {
                double ratio = (double) k / (double) K;
                calculatePhongModel(d, ray, i,
                                    obj, n, R * ratio, G * ratio, B * ratio, &tmp_col);
            }
        }
    }

    if (obj->alphaMap) alphaMap(obj->alphaMap, a, b, &obj->alpha);

    // no light source in the scene
    if (d.pw == -1) {
        return;
    }

    // ambient component
    obj->alb.ra = 0.3;
    tmp_col.R += obj->alb.ra * R;
    tmp_col.G += obj->alb.ra * G;
    tmp_col.B += obj->alb.ra * B;

    // Global component
    if (depth <= MAX_DEPTH) {

        // refraction

        if (obj->alpha + 1e-6 < 1) {
            tmp_col.R *= obj->alpha;
            tmp_col.G *= obj->alpha;
            tmp_col.B *= obj->alpha;

            double r_idx1 = ray->insideOut ? 1.0 : obj->r_index;
            double r_idx2 = ray->insideOut ? obj->r_index : 1.0;

            struct point3D *op_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);
            double cos_theta1 = dot(op_d, n);
            double sin_theta1 = sqrt(1 - pow(cos_theta1, 2));
            double sin_theta2 = (double) (r_idx1 / r_idx2) * sin_theta1;

            // has refraction , not total reflection
            if (sin_theta2 < 1 && sin_theta2 > 0) {
                double n21 = r_idx1 / r_idx2;
                double dot_product = -dot(n, &ray->d);
                double tmp = sqrt(1 - n21 * n21 * (1 - dot_product * dot_product));

                struct point3D *refract_d = newPoint(n21 * (dot_product * n->px + ray->d.px) - tmp * n->px,
                                                     n21 * (dot_product * n->py + ray->d.py) - tmp * n->py,
                                                     n21 * (dot_product * n->pz + ray->d.pz) - tmp * n->pz);
                normalize(refract_d);

                struct colourRGB refracted_color;
                struct ray3D *refracted_ray = newRay(p, refract_d);
                refracted_ray->insideOut = 1 - ray->insideOut;
                rayTrace(refracted_ray, depth + 1, &refracted_color, obj);

                tmp_col.R += (1 - obj->alpha) * R * refracted_color.R;
                tmp_col.G += (1 - obj->alpha) * B * refracted_color.G;
                tmp_col.B += (1 - obj->alpha) * G * refracted_color.B;

                free(refract_d);
                free(refracted_ray);
            }

            free(op_d);
        }
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
            initRay(&mirror_ray, p, &mirror_ray.d, 1);

            // Trace the ray and add its result
            rayTrace(&mirror_ray, depth + 1, &mirror_result, obj);
            tmp_col.R += mirror_result.R * obj->alb.rg;
            tmp_col.G += mirror_result.G * obj->alb.rg;
            tmp_col.B += mirror_result.B * obj->alb.rg;
        }
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

    // traverse all the object that is not itself or light source
    for (struct object3D *i = object_list; i != NULL; i = i->next) {
        if (i != Os && !i->isLightSource){
            struct point3D tmp_p, tmp_n;
            double tmp_a, tmp_b;

            //find intersection
            (i->intersect)(i, ray, &lambda_min, &tmp_p, &tmp_n, &tmp_a, &tmp_b);

            // replace first hit with the smallest positive lambda
            if (lambda_min > 0 && (lambda_min < *lambda || *lambda <= 0)) {
                *lambda = lambda_min;
                memcpy(p, &tmp_p, sizeof(struct point3D));
                memcpy(n, &tmp_n, sizeof(struct point3D));
                *a = tmp_a;
                *b = tmp_b;
                *obj = i;
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
        ray->insideOut = Os == obj ? 0 : 1;
        rtShade(obj, &p, &n, ray, depth, a, b, &I);
        // colRGB in  [0, 1]
        col->R = max(0, I.R);
        col->G = max(0, I.G);
        col->B = max(0, I.B);
    } else { // hit nothing
        col->R = 0;
        col->G = 0;
        col->B = 0;
    }
}

void calculatePixel(struct point3D *pc, struct colourRGB *col, struct view *cam, struct colourRGB *background) {
    // init variable
    struct point3D d;
    struct ray3D ray;
    pc->pz = cam->f;
    pc->pw = 1;

    // Convert pc into pw and calculate the direction of ray
    matVecMult(cam->C2W, pc); // pw = M_cw' * pc
    memcpy(&d, pc, sizeof(struct point3D));
    subVectors(&cam->e, &d); // d = pw - e
    normalize(&d);

    // Initialize the ray
    initRay(&ray, &cam->e, &d, 1);
    normalize(&ray.d);

    // Trace the ray of the pixel
    memcpy(col, background, sizeof(struct colourRGB));
    rayTrace(&ray, 1, col, NULL);
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
    int antialiasing_k = 5;
    struct colourRGB sample_RGB;
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

    // create the photon map image for each object except for the light source
    for (struct object3D *diff_obj = object_list; diff_obj != NULL; diff_obj = diff_obj->next) {
        if (!diff_obj->isLightSource && diff_obj->alb.rd != 0) { // is not a light source and is a diffuse object
            diff_obj->photonMap = newImage(sx, sx);
            diff_obj->photonMap->rgbdata = (double *)realloc(diff_obj->photonMap->rgbdata, sizeof(double) * sx * sx * 3);
        }
    }

    int num_rays = 100000;
    struct ray3D random_ray;
    for (struct object3D *ls_obj = object_list; ls_obj != NULL; ls_obj = ls_obj->next) {
        if (ls_obj->isLightSource) {
            for (int k = 0; k < num_rays; ++k) {
                ls_obj->initRandRay(ls_obj, &random_ray);
                forwardPassTrace(&random_ray, 1, ls_obj, ls_obj->col.R, ls_obj->col.G, ls_obj->col.B, sx);
            }
        }
    }

    fprintf(stderr, "Rendering row: ");
#pragma omp parallel for schedule(dynamic, 16) shared(antialiasing_k, rgbIm, object_list, cam, background) private(j)
    for (j = 0; j < sx; j++)        // For each of the pixels in the image
    {
        fprintf(stderr, "%d/%d, ", j, sx);
#pragma omp parallel for private(pc, d, ray, col, sample_RGB, i)
        for (i = 0; i < sx; i++) {
            ///////////////////////////////////////////////////////////////////
            // TO DO - complete the code that should be in this loop to do the
            //         raytracing!
            ///////////////////////////////////////////////////////////////////
            // initialize the final pixel color
            col.R = 0;
            col.G = 0;
            col.B = 0;

            if (antialiasing) {
                // 5-sampled anti-aliasing
                for (int k = 0; k < 5; ++k) {
                    // Calculate sample point in camera coordinate
                    // the sampling point should be with in [px-0.5, px+0.5) on x-coord
                    //                                      [py-0.5, py+0.5) on y-coord
                    pc.px = cam->wl + du * (i + ((double) rand()/ RAND_MAX) - 0.5);
                    pc.py = cam->wt + dv * (j + ((double) rand()/ RAND_MAX) - 0.5);

                    calculatePixel(&pc, &sample_RGB, cam, &background);
                    col.R += sample_RGB.R;
                    col.G += sample_RGB.G;
                    col.B += sample_RGB.B;
                }
                col.R = col.R / antialiasing_k;
                col.G = col.G / antialiasing_k;
                col.B = col.B / antialiasing_k;
            } else {
                // Calculate pc in camera coordinate
                pc.px = cam->wl + i * du;
                pc.py = cam->wt + j * dv;

                calculatePixel(&pc, &col, cam, &background);
            }

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

