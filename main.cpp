#include <iostream>
#include <parasolid.h>

#include <vector>
#include <map>
#include <sstream>
#include <string>
#include <iomanip>


bool is_psbody(int id) {
    PK_ERROR_t err = PK_ERROR_no_errors;
    PK_CLASS_t entity_class;
    err = PK_ENTITY_ask_class(id, &entity_class);
    return (err == 0) && (entity_class == PK_CLASS_body);
}

void verify_body(int body, std::stringstream& ss) {
    PK_ERROR_t err = PK_ERROR_no_errors;

    PK_BODY_check_o_t check_opts;
    PK_BODY_check_o_m(check_opts);
    int nfaults;
    PK_check_fault_t* faults;
    bool has_missing_geometry = false;
    bool has_invalid_state = false;
    bool has_corrupt_state = false;
    
    
    err = PK_BODY_check(body, &check_opts, &nfaults, &faults);
    if (nfaults) {

        for (int i = 0; i < nfaults; ++i) {
            auto check_state = faults[i].state;
            switch (check_state) {
            case PK_TOPOL_state_no_geom_c:
                has_missing_geometry = true;
                break;
            case PK_ENTITY_state_invalid_c:
                has_invalid_state = true;
                break;
            case PK_BODY_state_corrupt_c:
                has_corrupt_state = true;
                break;
            }
        }

        PK_MEMORY_free(faults);
    }


    PK_BODY_ask_topology_o_t options;
    PK_BODY_ask_topology_o_m(options);
    options.want_fins = PK_LOGICAL_false;
    PK_TOPOL_t* topols;
    PK_CLASS_t* classes;
    int n_topols;
    int n_relations;
    int* parents;
    int* children;
    PK_TOPOL_sense_t* senses;
    
    err = PK_BODY_ask_topology(body, &options, &n_topols, &topols, &classes, &n_relations, &parents, &children, &senses);

    PK_SURF_t surf;
    PK_CURVE_t curve;
    PK_POINT_t point;
    PK_CLASS_t entity_class;

    bool error_checking_topology = false;

    // missing geo counts
    int n_faces_no_geo = 0;
    int n_edges_no_geo = 0;
    int n_verts_no_geo = 0;
    // topo counts
    int n_regions = 0;
    int n_shells = 0;
    int n_faces = 0;
    int n_edges = 0;
    int n_loops = 0;
    int n_vertices = 0;

    // surf counts
    int nplane = 0;
    int ncyl = 0;
    int ncone = 0;
    int nsphere = 0;
    int ntorus = 0;
    int nbsurf = 0;
    int noffset = 0;
    int nfsurf = 0;
    int nswept = 0;
    int nspun = 0;
    int nblendsf = 0;
    
    // curve counts
    int n_line = 0;
    int n_circle = 0;
    int n_ellipse = 0;
    int n_bcurve = 0;
    int n_icurve = 0;
    int n_fcurve = 0;
    int n_spcurve = 0;
    int n_trcurve = 0;
    int n_cpcurve = 0;
    
    for (int i = 0; i < n_topols; ++i) {
        switch (classes[i]) {
        case PK_CLASS_part:
            break;
        case PK_CLASS_instance:
            break;
        case PK_CLASS_region:
            ++n_regions;
            break;
        case PK_CLASS_shell:
            ++n_shells;
            break;
        case PK_CLASS_face:
            ++n_faces;

            err = PK_FACE_ask_surf(topols[i], &surf);
            if (err != PK_ERROR_no_errors) {
                error_checking_topology = true;
                continue;
            }
            if (surf == PK_ENTITY_null) {
                ++n_faces_no_geo;
                continue;
            }
            err = PK_ENTITY_ask_class(surf, &entity_class);
            if (err != PK_ERROR_no_errors) {
                error_checking_topology = true;
                continue;
            }
            switch (entity_class) {
            case PK_CLASS_plane:
                ++nplane;
                break;
            case PK_CLASS_cyl:
                ++ncyl;
                break;
            case PK_CLASS_cone:
                ++ncone;
                break;
            case PK_CLASS_sphere:
                ++nsphere;
                break;
            case PK_CLASS_torus:
                ++ntorus;
                break;
            case PK_CLASS_bsurf:
                ++nbsurf;
                break;
            case PK_CLASS_offset:
                ++noffset;
                break;
            case PK_CLASS_fsurf:
                ++nfsurf;
                break;
            case PK_CLASS_swept:
                ++nswept;
                break;
            case PK_CLASS_spun:
                ++nspun;
                break;
            case PK_CLASS_blendsf:
                ++nblendsf;
                break;
            }

            break;
        case PK_CLASS_edge:
            ++n_edges;

            err = PK_EDGE_ask_curve(topols[i], &curve);
            if (err != PK_ERROR_no_errors) {
                error_checking_topology = true;
                continue;
            }
            if (curve == PK_ENTITY_null) {
                ++n_edges_no_geo;
                continue;
            }
            err = PK_ENTITY_ask_class(curve, &entity_class);
            if (err != PK_ERROR_no_errors) {
                error_checking_topology = true;
                continue;
            }
            switch (entity_class) {
            case PK_CLASS_line:
                ++n_line;
                break;
            case PK_CLASS_circle:
                ++n_circle;
                break;
            case PK_CLASS_ellipse:
                ++n_ellipse;
                break;
            case PK_CLASS_bcurve:
                ++n_bcurve;
                break;
            case PK_CLASS_icurve:
                ++n_icurve;
                break;
            case PK_CLASS_fcurve:
                ++n_fcurve;
                break;
            case PK_CLASS_spcurve:
                ++n_spcurve;
                break;
            case PK_CLASS_trcurve:
                ++n_trcurve;
                break;
            case PK_CLASS_cpcurve:
                ++n_cpcurve;
                break;
            }


            break;
        case PK_CLASS_fin:
            break;
        case PK_CLASS_loop:
            ++n_loops;
            break;
        case PK_CLASS_vertex:
            ++n_vertices;
            err = PK_VERTEX_ask_point(topols[i], &point);
            if (err != PK_ERROR_no_errors) {
                error_checking_topology = true;
                continue;
            }
            if (point == PK_ENTITY_null) {
                ++n_verts_no_geo;
                continue;
            }
            break;
        }
    }
    PK_TOPOL_eval_mass_props_o_t mp_opts;
    PK_TOPOL_eval_mass_props_o_m(mp_opts);
    double amount, mass, c_of_g[3], m_of_i[9], periphery;
    err = PK_TOPOL_eval_mass_props(1, &body, 0.999, &mp_opts, &amount, &mass, c_of_g, m_of_i, &periphery);
    
    bool error_computing_mass_properties = false;
    
    if (err != PK_ERROR_no_errors) {
        error_computing_mass_properties = true;
        amount = 0;
        mass = 0;
        for (int i = 0; i < 3; ++i) {
            c_of_g[i] = 0;
        }
        for (int i = 0; i < 9; ++i) {
            m_of_i[i] = 0;
        }
        periphery = 0;
    }

    // Bounding Box

    PK_BOX_t box;
    err = PK_TOPOL_find_box(body, &box);
    bool error_finding_bounding_box = false;
    if (err != PK_ERROR_no_errors) {
        error_finding_bounding_box = true;
        for (int i = 0; i < 6; ++i) {
            box.coord[i] = 0;
        }
    }

    // NA Bounding Box (may be better for distinguishing parts?)
    PK_TOPOL_find_nabox_o_t na_opts;
    PK_TOPOL_find_nabox_o_m(na_opts);
    na_opts.quality = PK_NABOX_quality_improved_c;
    PK_NABOX_sf_t nabox;
    bool error_finding_na_box = false;
    err = PK_TOPOL_find_nabox(1, &body, NULL, &na_opts, &nabox);
    if (err != PK_ERROR_no_errors) {
        error_finding_na_box = true;
        for (int i = 0; i < 3; ++i) {
            nabox.basis_set.axis.coord[i] = 0;
            nabox.basis_set.location.coord[i] = 0;
            nabox.basis_set.ref_direction.coord[i] = 0;
            nabox.box.coord[i] = 0;
            nabox.box.coord[i + 3] = 0;
        }
    }

    // |error_outputs| = 7
    std::vector<bool> error_outputs{ 
        has_corrupt_state, 
        has_invalid_state, 
        has_missing_geometry, 
        error_checking_topology, 
        error_finding_bounding_box, 
        error_finding_na_box, 
        error_computing_mass_properties
    };
    // |error_counts| = 4
    std::vector<int> error_counts{
        nfaults,
        n_faces_no_geo,
        n_edges_no_geo,
        n_verts_no_geo
    };

    // type_counts = {#nodes, #edges, topo counts, surf counts, curve counts}
    // size = 27
    std::vector<int> type_counts{
        n_topols,
        n_relations,
        n_regions,
        n_shells,
        n_faces,
        n_edges,
        n_loops,
        n_vertices,
        nplane,
        ncyl,
        ncone,
        nsphere,
        ntorus,
        nbsurf,
        noffset,
        nfsurf,
        nswept,
        nspun,
        nblendsf,
        n_line,
        n_circle,
        n_ellipse,
        n_bcurve,
        n_icurve,
        n_fcurve,
        n_spcurve,
        n_trcurve,
        n_cpcurve
    };

    // geom_props = {bb, nabb, mass_props}
    // size = 36
    std::vector<double> geom_props;
    for (int i = 0; i < 6; ++i) {
        geom_props.push_back(box.coord[i]);
    }
    for (int i = 0; i < 3; ++i) {
        geom_props.push_back(nabox.basis_set.axis.coord[i]);
    }
    for (int i = 0; i < 3; ++i) {
        geom_props.push_back(nabox.basis_set.location.coord[i]);
    }
    for (int i = 0; i < 3; ++i) {
        geom_props.push_back(nabox.basis_set.ref_direction.coord[i]);
    }
    for (int i = 0; i < 6; ++i) {
        geom_props.push_back(nabox.box.coord[i]);
    }
    geom_props.push_back(amount); 
    geom_props.push_back(mass);
    for (int i = 0; i < 3; ++i) {
        geom_props.push_back(c_of_g[i]);
    }
    for (int i = 0; i < 9; ++i) {
        geom_props.push_back(m_of_i[i]);
    }
    geom_props.push_back(periphery);
    
    // overall output length = 38 int + 36 doubles = 74
    
    for (auto error : error_outputs) {
        ss << "," << ((int)error);
    }
    for (auto count : error_counts) {
        ss << "," << count;
    }
    for (auto count : type_counts) {
        ss << "," << count;
    }
    // Maximize output precision for double properties
    ss << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    for (auto prop : geom_props) {
        ss << "," << prop;
    }

}

std::string fingerprint_part(char* path) {
    ensure_parasolid_session();
    PK_PART_receive_o_t receive_opts;
    PK_PART_receive_o_m(receive_opts);
    receive_opts.transmit_format = PK_transmit_format_text_c;
    int n_parts = 0;
    PK_PART_t* parts = NULL;
    PK_ERROR_t err = PK_ERROR_no_errors;
    err = PK_PART_receive(path, &receive_opts, &n_parts, &parts);

    std::vector<int> ps_bodies;
    bool can_read_file = false;

    if (err == 0) {
        can_read_file = true;
        for (int i = 0; i < n_parts; ++i) {
            if (is_psbody(parts[i])) {
                ps_bodies.push_back(parts[i]);
            }
        }
    }
    else {
        n_parts = 0;
    }

    int num_bodies = ps_bodies.size();

    std::stringstream ss;
    ss << int(can_read_file) << "," << n_parts << "," << num_bodies;

    if (num_bodies == 1) {
        verify_body(ps_bodies[0], ss);
    }
    else { // Pad with 0s if part is invalid for number of bodies
        for (int i = 0; i < 74; ++i) {
            ss << "," << 0;
        }
    }


    PK_MEMORY_free(parts);

    return ss.str();
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: verisolid path1 [path2 ...]" << std::endl;
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        std::string part_info = fingerprint_part(argv[i]);
        std::cout << part_info << std::endl;
    }
   
    return 0;
}
