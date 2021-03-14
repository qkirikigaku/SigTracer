#include "SigTracer.h"

void run_VB(string data_directory, string sample_name, int start_J,
        int iteration, vector<string> ref_sig, string ABS_path, string mode, bool with_noise){
    while(true){
        SigTracer st(data_directory, sample_name, start_J, iteration, ref_sig, ABS_path, mode, with_noise);
        st.load_data();
        st.run();
        if(st.exit_call == false){
            st.write_result();
            break;
        }
    }
}
