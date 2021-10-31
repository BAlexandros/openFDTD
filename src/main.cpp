#include "../include/grid1D.hpp"
#include "../include/imgui.h"
#include "../include/imgui_impl_opengl3.h"
#include "../include/imgui_impl_glfw.h"
#include <GLFW/glfw3.h>
#include <omp.h>

#ifndef __has_include
static_assert(false, "__has_include not supported");
#else
#  if __cplusplus >= 201703L && __has_include(<filesystem>)
#    include <filesystem>
namespace fs = std::filesystem;
#  elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#  elif __has_include(<boost/filesystem.hpp>)
#    include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#  endif
#endif

static void glfw_error_callback(int error, const char* description);

int main(void)
{

  // Setup window
  glfwSetErrorCallback(glfw_error_callback);
  if (!glfwInit())
    return 1;
  const char* glsl_version = "#version 130";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

  // Create window with graphics context
  GLFWwindow* window = glfwCreateWindow(1280, 720, "openFDTD", NULL, NULL);
  if (window == NULL)
    return 1;
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1); // Enable vsync

  // Setup Dear ImGui context
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO(); (void)io;
  ImFontConfig fontconfig;
  fontconfig.OversampleH = 6;
  fontconfig.OversampleV = 6;
  fontconfig.GlyphExtraSpacing.x = 1.0f;
  io.Fonts->AddFontDefault();
  ImFont* font1 = io.Fonts->AddFontFromFileTTF("fonts/source-sans-pro/SourceSansPro-Semibold.otf", 15, &fontconfig);
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;       // Enable Keyboard Controls

  // Setup Dear ImGui style
  // ImGui::StyleColorsDark();
  ImGui::StyleColorsClassic();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  // Our state
  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
  
  static bool cleanup_temp_frames_at_exit = true;

  // Main loop
  while (!glfwWindowShouldClose(window))
  {
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    
    ImGui::PushFont(font1);

    ImGui::Begin("FDTD Settings",NULL,ImGuiWindowFlags_MenuBar);

    // Menu Bar
    static bool show_exit_confirm = false;
    ImGui::BeginMenuBar();
    if (ImGui::MenuItem("Exit")){
      show_exit_confirm = true;
    }
    ImGui::EndMenuBar();

    // Exit dialog confirmation
    if (show_exit_confirm){
      ImGui::OpenPopup("Exit?");

      ImVec2 center = ImGui::GetMainViewport()->GetCenter();
      ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

      if (ImGui::BeginPopupModal("Exit?", NULL, ImGuiWindowFlags_AlwaysAutoResize))
      {
        ImGui::Text("Are you sure you want to exit?\n\n");
        ImGui::Separator();
        if (ImGui::Button("Yes", ImVec2(120, 0))) { ImGui::CloseCurrentPopup(); glfwSetWindowShouldClose(window, 1); }
        ImGui::SameLine();
        if (ImGui::Button("No", ImVec2(120, 0))) { show_exit_confirm = false; ImGui::CloseCurrentPopup(); }
        ImGui::SetItemDefaultFocus();
        ImGui::EndPopup();
      }
    }


    /*****************************
     **   DIMENSION SELECTION   **
     *****************************/
    ImGui::Text("Dimensions");
    static int selected_dimN = 1;
    ImGui::RadioButton("1D", &selected_dimN, 1); ImGui::SameLine();
    ImGui::RadioButton("2D", &selected_dimN, 2); ImGui::SameLine();
    ImGui::RadioButton("3D", &selected_dimN, 3);
    static Grid1D *g1 = new Grid1D;

    /*****************************
     **   MULTITHREADING        **
     *****************************/

    ImGui::Separator();
    static bool parallelism_enabled_checkbox = false;
    static int max_num_threads = omp_get_max_threads();
    static int enabled_thread_num = 1;

    ImGui::Text("Multithreading");
    ImGui::Checkbox("Enable OpenMP", &parallelism_enabled_checkbox);
    if (parallelism_enabled_checkbox){
      ImGui::SameLine();
      ImGui::SetNextItemWidth(100);
      if (ImGui::SliderInt("Thread number", &enabled_thread_num, 1, max_num_threads, "%d")){
        omp_set_num_threads(enabled_thread_num);
      }
    }

    ImGui::Separator();

    ImGui::BeginGroup();
    ImGui::BeginGroup();

    /*************************
     **   SOURCE SETTINGS   **
     *************************/
    ImGui::Text("Source Settings");

    const char* sources_list[3] = {"Gaussian", "Sinusoidal", "Ricker"};
    static int selected_source = 0;
    static double selected_E0 = 1e0;
    static double selected_f  = 1.e6;
    if (selected_dimN == 1){
      ImGui::SetNextItemWidth(100);
      ImGui::Combo("Source Type", &selected_source, sources_list, 3);
    }

    ImGui::SetNextItemWidth(100);
    ImGui::InputDouble("Source Amplitude", &selected_E0, 0.0f, 0.0f, "%g");

    ImGui::SetNextItemWidth(100);
    ImGui::InputDouble("Source Frequency", &selected_f, 0.0f, 0.0f, "%.2e");

    /*************************
     **   GRID SETTINGS     **
     *************************/
    static int selected_Nx1D = 200;
    static int selected_Nt1D = 500;
    static int selected_Nl1D = 100;
    static double selected_S1D  = 1.0;
    ImGui::Text("Grid Settings");
    if (selected_dimN == 1){

      ImGui::PushItemWidth(50);
      ImGui::InputInt("Number of spatial steps in grid", &selected_Nx1D,0,0);
      if (selected_Nx1D <= 0 ) { selected_Nx1D = 1; };

      ImGui::InputInt("Number of time steps to simulate", &selected_Nt1D,0,0);
      if (selected_Nt1D <= 0 ) { selected_Nt1D = 1; };

      ImGui::InputInt("Number of spatial steps per wavelength", &selected_Nl1D,0,0);
      if (selected_Nl1D <= 0 ) { selected_Nl1D = 1; };

      ImGui::InputDouble("Courant number", &selected_S1D, 0.0f, 0.0f, "%.2f");
      if (selected_S1D <= 0 ) { selected_S1D = 1; };
      if (selected_S1D > 1) { ImGui::SameLine(); ImGui::TextColored(ImVec4(1,1,0,1),"WARNING: Courant number too large\nSimulation unstable");}
      ImGui::PopItemWidth();
    }


    /*****************************
     **   MATERIAL SETTINGS     **
     *****************************/

    ImGui::Text("Materials");
    static bool show_material_form = false;
    static std::vector<GridMat> selectedMats1D;
    static int selected_bounds[2] {0,0};
    static int selected_mat = 0;

    // Get the list of material names available in the database
    static char** mat_input_list = new char*[materialdb.size()];
    {
    int i = 0;
    for (auto const& x : materialdb)
    {
      mat_input_list[i] = new char[64];
      strcpy(mat_input_list[i], x.first.c_str());
      i++;
    }
    }

    if (selected_dimN == 1){
      if (ImGui::Button("Add material")){
        show_material_form = show_material_form ? false : true;
      }
    }
    if (show_material_form){
      if (selected_dimN == 1){
        ImGui::PushItemWidth(150);
        ImGui::InputInt2("Material Bounds", selected_bounds);
        bool bounds_valid = selected_bounds[0] > 0 &&
                            selected_bounds[1] < selected_Nx1D &&
                            selected_bounds[0] < selected_bounds[1];
        if (!bounds_valid){
          ImGui::SameLine(); 
          ImGui::TextColored(ImVec4(1,1,0,1),"!!");
          if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Boundaries are not valid");
        }
        ImGui::Combo("Material Type", &selected_mat, mat_input_list, materialdb.size());
        if (ImGui::Button("Add")){
          if (bounds_valid){
            g1->add_material(selected_bounds[0],selected_bounds[1],mat_input_list[selected_mat]);
            selected_bounds[0] = 0;
            selected_bounds[1] = 0;
            selected_mat = 0;
          }
        }
        ImGui::PopItemWidth();
      }
    }

    ImGui::EndGroup();
    ImGui::SameLine();
    ImGui::BeginGroup();

    /*****************************
     **   OUTPUT FILES          **
     *****************************/
    
    ImGui::Text("Output files");
    static char selected_field_data_fname[64]       = "data/field.dat"; 
    static char selected_spect_data_fname[64]       = "data/spectrum.dat"; 
    static char selected_field_animation_fname[64]  = "gallery/field.mp4"; 
    static char selected_spect_output_fname[64]     = "gallery/spectrum.png"; 
    static int  selected_animation_start_index = 0;
    static int  selected_animation_end_index = selected_Nt1D;
    double nyquist_max_f = 0.5/(selected_S1D/selected_f/selected_Nl1D);
    static double  selected_spectrum_max = nyquist_max_f;
    ImGui::PushItemWidth(150);
    ImGui::InputText("Field data filename",       selected_field_data_fname,      64, ImGuiInputTextFlags_CharsNoBlank);
    if (fs::exists(selected_field_data_fname)){
      ImGui::SameLine(); 
      ImGui::TextColored(ImVec4(1,1,0,1),"!!");
      if (ImGui::IsItemHovered())
        ImGui::SetTooltip("File name already exists\nIt will be overwritten");
    }
    ImGui::InputText("Spectrum data filename",    selected_spect_data_fname,      64, ImGuiInputTextFlags_CharsNoBlank);
    if (fs::exists(selected_spect_data_fname)){
      ImGui::SameLine(); 
      ImGui::TextColored(ImVec4(1,1,0,1),"!!");
      if (ImGui::IsItemHovered())
        ImGui::SetTooltip("File name already exists\nIt will be overwritten");
    }
    ImGui::InputText("Field animation filename",  selected_field_animation_fname, 64, ImGuiInputTextFlags_CharsNoBlank);
    if (fs::exists(selected_field_animation_fname)){
      ImGui::SameLine(); 
      ImGui::TextColored(ImVec4(1,1,0,1),"!!");
      if (ImGui::IsItemHovered())
        ImGui::SetTooltip("File name already exists\nIt will be overwritten");
    }
    ImGui::InputText("Spectrum output filename",  selected_spect_output_fname, 64, ImGuiInputTextFlags_CharsNoBlank);
    if (fs::exists(selected_spect_output_fname)){
      ImGui::SameLine(); 
      ImGui::TextColored(ImVec4(1,1,0,1),"!!");
      if (ImGui::IsItemHovered())
        ImGui::SetTooltip("File name already exists\nIt will be overwritten");
    }
    ImGui::PopItemWidth();
    ImGui::PushItemWidth(50);
    ImGui::InputInt("Animation first index", &selected_animation_start_index,0,0);
    ImGui::SameLine();
    if (ImGui::SmallButton("Make start")){
      selected_animation_start_index = 0;
    }
    ImGui::InputInt("Animation last  index", &selected_animation_end_index,0,0);
    ImGui::SameLine();
    if (ImGui::SmallButton("Make end")){
      selected_animation_end_index = selected_Nt1D;
    }
    selected_animation_start_index = (selected_animation_start_index > 0 ||
                                      selected_animation_start_index < selected_Nt1D ||
                                      selected_animation_start_index < selected_animation_end_index) 
                                      ? selected_animation_start_index : 0;
    selected_animation_end_index = (selected_animation_end_index > 0 ||
                                      selected_animation_end_index < selected_Nt1D ||
                                      selected_animation_start_index < selected_animation_end_index) 
                                      ? selected_animation_end_index : selected_Nt1D;
    ImGui::InputDouble("Maximum spectrum f (Hz)",&selected_spectrum_max,0,0,"%g");
    selected_spectrum_max = (selected_spectrum_max <= 0 || selected_spectrum_max > nyquist_max_f) ? nyquist_max_f : selected_spectrum_max;
    ImGui::SameLine();
    if (ImGui::SmallButton("Make max")){
      selected_spectrum_max = nyquist_max_f;
    }
    if (ImGui::IsItemHovered()){
      ImGui::SetTooltip("Nyquist sampling maximum\nf = %g Hz",nyquist_max_f);
    }
    ImGui::PopItemWidth();
    ImGui::Checkbox("Delete frames on exit", &cleanup_temp_frames_at_exit);


    /*****************************
     **   SIMULATION EXECUTION  **
     *****************************/

    // If the setup succeeds, this flag indicates that 
    // the simulation can be run
    static bool grid_ready = false;
    static bool data_generated = false;
    // Setup the grid according to the user input
    if (ImGui::Button("Setup Simulation")){
      if ( parallelism_enabled_checkbox){
        g1->enable_parallelism();
        g1->setThreadNumber(enabled_thread_num);
      } else {
        g1->disable_parallelism();
      }

      g1->setE0(selected_E0);
      g1->setFreq(selected_f);
      g1->setNx(selected_Nx1D);
      g1->setNl(selected_Nl1D);
      g1->setSc(selected_S1D);
      g1->setNt(selected_Nt1D);
      if (strcmp(sources_list[selected_source], "Gaussian") == 0){
        g1->setSourceType(Source_t::Gaussian);
      } else if (strcmp(sources_list[selected_source], "Sinusoidal") == 0) {
        g1->setSourceType(Source_t::Sinusoidal);
      } else if (strcmp(sources_list[selected_source], "Ricker") == 0) {
        g1->setSourceType(Source_t::Ricker);
      }
      g1->setFieldProgFname(selected_field_data_fname);
      g1->setSpectrumFname(selected_spect_data_fname);
      g1->setSpectOutFname(selected_spect_output_fname);
      for (size_t i = 0; i < g1->material_list.size(); i++){
        if (g1->material_list[i].x2 >= selected_Nx1D){
          g1->material_list[i].x2 = selected_Nx1D - 1;
        }
      }
      g1->animation_start_index = selected_animation_start_index;
      g1->animation_end_index   = selected_animation_end_index - 1;
      g1->spectrum_max_index = selected_spectrum_max;
      grid_ready = true;
    }
    if (grid_ready ){
      ImGui::SameLine();
      if (ImGui::Button("Run Simulation")){
        fs::create_directory("temp_openFDTD_frames");
        g1->run_simulation();

        char ffmpeg_command[512];
        sprintf(ffmpeg_command,\
                    "ffmpeg -y -hide_banner -loglevel error "
                    "-r 60 -f image2 -s 1920x1080 "
                    " -start_number %d "
                    "-i temp_openFDTD_frames/frame-%%05d.png "
                    "-vframes %d "
                    "-vcodec libx264 -crf 25 -pix_fmt yuv420p %s", \
               selected_animation_start_index, 
               selected_animation_end_index - selected_animation_start_index,\
               selected_field_animation_fname);
        system(ffmpeg_command);
        data_generated = true;
      }
    }

    // The grid can be reset after it is created or the simulation run
    if (grid_ready){
      if (ImGui::Button("Reset grid")){
        delete g1;
        g1 = new Grid1D();
        grid_ready = false;
        data_generated = false;
      }
    }
    
    if (data_generated){
      ImGui::SameLine();
      if (ImGui::Button("Remake spectrum plot")){
        g1->spectrum_max_index = selected_spectrum_max;
        g1->makeSpectrumPlot(); 
      }
    }

    ImGui::EndGroup();
    ImGui::EndGroup();

    /*****************************
     **   Material List Review  **
     *****************************/

    if (!show_material_form){
      ImVec2 cursor = ImGui::GetCursorPos();
      ImGui::SetCursorPos(ImVec2(cursor[0], cursor[1] + 75));

    }

    static bool show_mats_in_grid = true;
    ImGui::Text("Materials in grid");
    ImGui::SameLine();
    if (show_mats_in_grid){
      if (ImGui::SmallButton("Hide")){
        show_mats_in_grid = show_mats_in_grid ? false : true;
      }
    } else {
      if (ImGui::SmallButton("Show")){
        show_mats_in_grid = show_mats_in_grid ? false : true;
      }
    }
    ImGui::Separator();


    // Generate a table of the materials in the grid
    if (show_mats_in_grid) {
      if (selected_dimN == 1){
        static ImGuiTableFlags mat_table_flags = \
                                                 ImGuiTableFlags_BordersV |          \
                                                 ImGuiTableFlags_BordersOuterH |     \
                                                 ImGuiTableFlags_RowBg |             \
                                                 ImGuiTableFlags_ContextMenuInBody | \
                                                 ImGuiTableFlags_NoHostExtendX |     \
                                                 ImGuiTableFlags_SizingFixedFit;
        ImGui::BeginTable("Material List",5, mat_table_flags);
        ImGui::TableSetupColumn("ID");
        ImGui::TableSetupColumn("Left Boundary");
        ImGui::TableSetupColumn("Right Boundary");
        ImGui::TableSetupColumn("Material");
        ImGui::TableHeadersRow();

        for (unsigned int i = 0; i < g1->material_list.size(); i++){
          ImGui::TableNextColumn();
          ImGui::Text("%3u",i);
          ImGui::TableNextColumn();
          ImGui::Text("%d",g1->material_list[i].x1);
          ImGui::TableNextColumn();
          ImGui::Text("%d",g1->material_list[i].x2);
          ImGui::TableNextColumn();
          ImGui::Text("%s",g1->material_list[i].matname.c_str());
          ImGui::TableNextColumn();
          if (ImGui::Button("Delete")){
            g1->material_list.erase(g1->material_list.begin()+i);
          }
        }
        ImGui::EndTable();
      }
    }
    ImGui::PopFont();

    ImGui::End();

    // Rendering
    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(window);
  }
  
  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  glfwDestroyWindow(window);
  glfwTerminate();

  if (cleanup_temp_frames_at_exit){
    fs::remove_all("./temp_openFDTD_frames");
  }

  return 0;
}

static void glfw_error_callback(int error, const char* description)
{
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}
