// cromq_wrapper.cpp

#include "cromq_wrapper.hpp"


CROMq_wrapper::CROMq_wrapper(std::string name_in,
                             int num_x_in,
                             int xdim_in,
                             double rd_param_in,
                             double R_enc_in,
                             int r_dec_num_in,
                             std::vector<double> R_dec_array_in,
                             std::string svd_folder_in,
                             std::string log_folder_in,
                             bool verbose_in) {
    name = name_in;
    num_x = num_x_in;
    xdim = xdim_in;
    num_line = 0;
    
    rd_param = rd_param_in;
    R_enc = R_enc_in;
    
    verbose = verbose_in;
    r_dec_num = r_dec_num_in;

    R_enc_per_col.resize(num_x);
    R_dec_array.resize(r_dec_num);

    for (int dec_idx=0; dec_idx<r_dec_num; dec_idx++) {
        R_dec_array[dec_idx] = R_dec_array_in[dec_idx];
    }
    mu.resize(num_x);
    std_array.resize(num_x);
    v_mat.resize(num_x, std::vector<double>(num_x));
    cov.resize(num_x, std::vector<double>(num_x));
    subfnames.resize(num_x);

    for (int row_idx=0; row_idx<num_x; row_idx++) {
        mu[row_idx] = 0;
        for (int col_idx=0; col_idx<num_x; col_idx++) {
            cov[row_idx][col_idx] = 0;
        }
    }

    svd_folder = svd_folder_in;
    log_folder = log_folder_in;

    // Create log folder and log file
    std::string log_fname = get_log_fname(log_folder, R_enc, rd_param, xdim);

    log_file.open(log_fname);
    log_file << log_fname << std::endl;

    // set svd parameter filenames
    int status = mkdir(svd_folder.c_str(), 0777);
    std::string mu_fname = svd_folder + "/mu_" + name + ".txt";
    std::string v_mat_fname = svd_folder + "/v_mat_" + name + ".txt";
    std::string std_array_fname = svd_folder + "/std_array_" + name + ".txt";
}

CROMq_wrapper::~CROMq_wrapper(){}

std::string CROMq_wrapper::get_subfname(int file_idx) {
    std::stringstream subfname_tmp;
    subfname_tmp.str(std::string());
    subfname_tmp << std::setw(4) << std::setfill('0') << file_idx;

    std::string subfname = "split_" + name + "_" + subfname_tmp.str() + ".subqscores";
    return subfname;
}

void CROMq_wrapper::run_svd() {
    int row_idx, col_idx;
    MatrixXf cov_eigen(num_x, num_x);
    MatrixXf sValues, vMatrix;

    for (row_idx=0; row_idx<num_x; row_idx++) {
        for (col_idx=0; col_idx<num_x; col_idx++) {
            cov_eigen(row_idx, col_idx) = cov[row_idx][col_idx];
        }
    }
    JacobiSVD<MatrixXf> svd(cov_eigen, ComputeThinU | ComputeThinV);
    sValues = svd.singularValues();
    vMatrix = svd.matrixV();
    for (row_idx=0; row_idx<num_x; row_idx++) {
        std_array[row_idx] = sqrt(sValues(row_idx));
        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat[row_idx][col_idx] = vMatrix(row_idx, col_idx);
        }
    }
}

void CROMq_wrapper::write_svd_params() {
    std::ofstream mu_file(mu_fname);
    std::ofstream v_mat_file(v_mat_fname);
    std::ofstream std_array_file(std_array_fname);

    int row_idx, col_idx;
    for (row_idx=0; row_idx<num_x; row_idx++) {
        mu_file << mu[row_idx] << std::endl;
        std_array_file << std_array[row_idx] << std::endl;
        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat_file << v_mat[row_idx][col_idx] << " ";
        }
        v_mat_file << std::endl;
    }
    mu_file.close();
    v_mat_file.close();
    std_array_file.close();
}

void CROMq_wrapper::read_svd_params() {
    int row_idx, col_idx;

    std::ifstream mu_file(mu_fname);
    std::ifstream std_array_file(std_array_fname);
    std::ifstream v_mat_file(v_mat_fname);

    if (verbose)
        std::cout << "Reading svd params from " << svd_folder << std::endl;

    for (row_idx=0; row_idx<num_x; row_idx++) {
        mu_file >> mu[row_idx];
    }
    mu_file.close();

    for (row_idx=0; row_idx<num_x; row_idx++) {
        std_array_file >> std_array[row_idx];
    }
    std_array_file.close();

    for (row_idx=0; row_idx<num_x; row_idx++) {
        for (col_idx=0; col_idx<num_x; col_idx++) {
            v_mat_file >> v_mat[row_idx][col_idx];
        }
    }
    v_mat_file.close();
}


void CROMq_wrapper::split_qscores(bool skip) {
    int line_idx = 0;
    int file_idx = 0;

    std::string fname = name + ".qscore";

    std::ifstream qscore_file(fname);
    std::string line;

    bool end_indicator = true;

    // split qscore file and get mu
    std::string subfname;
    while (end_indicator) {
        subfname = get_subfname(file_idx);
        subfnames[file_idx] = subfname;
        std::ofstream qscore_subfile(subfname);

        end_indicator = false;
        while(std::getline(qscore_file, line)) {
            qscore_subfile << line << std::endl;
            line_idx++;
            num_line++;
            for (int row_idx=0; row_idx<num_x; row_idx++) {
                mu[row_idx] += line.at(row_idx);
            }
            if (line_idx == xdim) {
                end_indicator = true;
                break;
            }
        }
        qscore_subfile.close();
        line_idx = 0;
        file_idx++;
    }
    num_file = file_idx-1;
    // XXX TBD take care of remainder (also when running cromq in the below)

    // normalize mu
    for (int row_idx=0; row_idx<num_x; row_idx++) {
        mu[row_idx] /= (double)num_line;
    }

    // reading from the beginning to compute the covariance matrix
    qscore_file.clear();
    qscore_file.seekg(0, std::ios::beg);

    std::vector<double> q_temp(num_x);

    while (std::getline(qscore_file, line)) {
        for (int row_idx=0; row_idx<num_x; row_idx++)
            q_temp[row_idx] = (line.at(row_idx) - mu[row_idx]);

        for (int row_idx=0; row_idx<num_x; row_idx++) {
            for (int col_idx=0; col_idx<=row_idx; col_idx++) {
                cov[row_idx][col_idx] += (q_temp[row_idx] * q_temp[col_idx]);
            }
        }
    }
    qscore_file.close();

    // normalize covariance matrix
    for (int row_idx=0; row_idx<num_x; row_idx++) {
        cov[row_idx][row_idx] /= (double)xdim;
        for (int col_idx=0; col_idx<row_idx; col_idx++) {
            cov[row_idx][col_idx] /= (double)xdim;
            cov[col_idx][row_idx] = cov[row_idx][col_idx];
        }
    }

    // run svd
    if (verbose)
        std::cout << "Running SVD" << std::endl;

    run_svd();

    // Write SVD parameters
    write_svd_params();
}

void CROMq_wrapper::run_encoders() {
    double enc_runtime = 0;
    std::clock_t run_time;

    // Run CROMq encoder
    for (int cromq_idx=0; cromq_idx<num_file; cromq_idx++) {
        run_time = std::clock();
        std::cout << "Processing " << subfnames[cromq_idx] << std::endl;
        CROMq_encoder enc(name,
                          subfnames[cromq_idx],
                          cromq_idx,
                          num_x,
                          xdim,
                          rd_param,
                          R_enc_per_col,
                          mu,
                          std_array,
                          v_mat,
                          verbose);
        enc.run();
        run_time = std::clock() - run_time;
        double runtime_sec = ((double)run_time / (double)CLOCKS_PER_SEC);
        std::cout << "Enc took " << runtime_sec << " seconds" << std::endl;
        log_file << "Enc took " << runtime_sec << " seconds" << std::endl;
        enc_runtime += runtime_sec;
    }
    std::cout << "Total: Enc took " << enc_runtime << " seconds" << std::endl << std::endl;
    log_file << "Total: Enc took " << enc_runtime << " seconds" << std::endl << std::endl;
}

void CROMq_wrapper::run_decoders() {
    double dec_runtime;
    std::clock_t run_time;
    double R_dec;
    std::vector<double> R_dec_per_col(num_x);

    // Run CROMq decoder
    for (int r_dec_idx=0; r_dec_idx<r_dec_num; r_dec_idx++) {
        dec_runtime = 0;
        run_time = std::clock();
        R_dec = R_dec_array[r_dec_idx];
        allocate_rate(std_array, R_dec_per_col, num_x, rd_param, R_dec);
        for (int cromq_idx=0; cromq_idx<num_file; cromq_idx++) {
            run_time = std::clock();
            std::cout << "Processing " << subfnames[cromq_idx];
            std::cout << " at rate = " << R_dec << std::endl;
            log_file << "Processing " << subfnames[cromq_idx];
            log_file << " at rate = " << R_dec << std::endl;
            CROMq_decoder dec(name,
                              subfnames[cromq_idx],
                              cromq_idx,
                              num_x,
                              xdim,
                              rd_param,
                              R_dec,
                              R_dec_per_col,
                              mu,
                              std_array,
                              v_mat,
                              verbose);
            dec.run();
            run_time = std::clock() - run_time;
            double runtime_sec = ((double)run_time / (double)CLOCKS_PER_SEC);
            std::cout << "Dec took " << runtime_sec << " seconds" << std::endl;
            log_file << "Dec took " << runtime_sec << " seconds" << std::endl;
            dec_runtime += runtime_sec;
        }
        std::cout << "Total: Dec took " << dec_runtime << " seconds" << std::endl;
        std::cout << std::endl << std::endl;
        log_file << "Total: Dec took " << dec_runtime << " seconds" << std::endl;
        log_file << std::endl << std::endl;
    }
}

void CROMq_wrapper::compute_distortions() {
    double distortion_avg;
    double R_dec;
    for (int dec_idx=0; dec_idx<r_dec_num; dec_idx++) {
        distortion_avg = 0;
        R_dec = R_dec_array[dec_idx];
        for (int cromq_idx=0; cromq_idx<num_file; cromq_idx++) {
            double distortion;
            std::string ofname = get_ofname(name, cromq_idx, R_dec);
            std::string ifname = get_subfname(cromq_idx);
            distortion = compute_distortion(ifname, ofname, num_x, xdim);
            std::cout << ofname << std::endl;
            std::cout << ifname << std::endl;
            std::cout << "Distortion (id " << cromq_idx << "): " << distortion << std::endl;
            log_file << "Distortion (id " << cromq_idx << "): " << distortion << std::endl;
            distortion_avg += distortion;
        }
        distortion_avg /= (double)(num_file);

        std::cout << "Average distortion : " << distortion_avg << std::endl;
        std::cout << std::endl << std::endl;
        log_file << "Average distortion : " << distortion_avg << std::endl;
        log_file << std::endl << std::endl;
    }
}

void CROMq_wrapper::run(bool encode, bool decode) {
    bool skip = true;
    if (verbose)
        std::cout << "splitting qscore file" << std::endl;
    split_qscores(skip);
    std::cout << "num_file = " << num_file << std::endl;
    log_file << "num_file = " << num_file << std::endl;
    allocate_rate(std_array, R_enc_per_col, num_x, rd_param, R_enc);

    if (encode)
        if (verbose)
            std::cout << "encoding" << std::endl;
        run_encoders();

    if (decode) {
        if (verbose)
            std::cout << "decoding" << std::endl;
        read_svd_params();
        run_decoders();
    }
    compute_distortions();
}
