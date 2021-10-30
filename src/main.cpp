#include "../include/grid1D.hpp"
#include "../include/imgui.h"
#include "../include/imgui_impl_opengl3.h"
#include "../include/imgui_impl_glfw.h"
#include <GLFW/glfw3.h>

static void glfw_error_callback(int error, const char* description);

int main(void)
{

  // Setup window
  glfwSetErrorCallback(glfw_error_callback);
  if (!glfwInit())
    return 1;
  const char* glsl_version = "#version 130";
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

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
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;       // Enable Keyboard Controls

  // Setup Dear ImGui style
  // ImGui::StyleColorsDark();
  ImGui::StyleColorsClassic();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);

  // Our state
  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

  // Main loop
  while (!glfwWindowShouldClose(window))
  {
    glfwPollEvents();

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    ImGui::Begin("FDTD Settings");

    /*****************************
     **   DIMENSION SELECTION   **
     *****************************/
    ImGui::Text("Dimensions");
    static int selected_dimN = 1;
    ImGui::RadioButton("1D", &selected_dimN, 1); ImGui::SameLine();
    ImGui::RadioButton("2D", &selected_dimN, 2); ImGui::SameLine();
    ImGui::RadioButton("3D", &selected_dimN, 3);
    static Grid1D *g1 = new Grid1D;

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

      ImGui::SetNextItemWidth(100);
      ImGui::InputInt("Number of spatial steps in grid", &selected_Nx1D, 1, 1);
      if (selected_Nx1D <= 0 ) { selected_Nx1D = 1; };

      ImGui::SetNextItemWidth(100);
      ImGui::InputInt("Number of time steps to simulate", &selected_Nt1D, 1, 1);
      if (selected_Nt1D <= 0 ) { selected_Nt1D = 1; };

      ImGui::SetNextItemWidth(100);
      ImGui::InputInt("Number of spatial steps per wavelength", &selected_Nl1D, 1, 1);
      if (selected_Nl1D <= 0 ) { selected_Nl1D = 1; };

      ImGui::SetNextItemWidth(100);
      ImGui::InputDouble("Courant number", &selected_S1D, 0.0f, 0.0f, "%.2f");
      if (selected_S1D <= 0 ) { selected_S1D = 1; };
      if (selected_S1D > 1) { ImGui::SameLine(); ImGui::TextColored(ImVec4(1,1,0,1),"WARNING: Courant number too large\nSimulation unstable");}
    }


    /*****************************
     **   MATERIAL SETTINGS     **
     *****************************/

    ImGui::Text("Materials");
    static bool show_material_form = false;
    static std::vector<GridMat> selectedMats1D;
    static int selected_bounds[2] {0,0};
    static int selected_mat = 0;

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
        ImGui::SetNextItemWidth(100);
        ImGui::InputInt2("Material Bounds", selected_bounds);
        ImGui::SetNextItemWidth(100);
        ImGui::Combo("Material Type", &selected_mat, mat_input_list, materialdb.size());
        if (ImGui::Button("Add")){
          if (selected_bounds[0] >= 0 &&
              selected_bounds[1] < selected_Nx1D &&
              selected_bounds[0] < selected_bounds[1]){
            g1->add_material(selected_bounds[0],selected_bounds[1],mat_input_list[selected_mat]);
            selected_bounds[0] = 0;
            selected_bounds[1] = 0;
            selected_mat = 0;
          }
        }
      }
    }


    ImGui::EndGroup();
    ImGui::SameLine();
    ImGui::BeginGroup();

    /*****************************
     **   OUTPUT FILES          **
     *****************************/
    
    ImGui::Text("Output files");
    static char selected_field_data_fname[64]       = "field.dat"; 
    static char selected_spect_data_fname[64]       = "spectrum.dat"; 
    static char selected_field_animation_fname[64]  = "field.gif"; 
    ImGui::PushItemWidth(100);
    ImGui::InputText("Field data filename",       selected_field_data_fname,      64, ImGuiInputTextFlags_CharsNoBlank);
    ImGui::InputText("Spectrum data filename",    selected_spect_data_fname,      64, ImGuiInputTextFlags_CharsNoBlank);
    ImGui::InputText("Field animation filename",  selected_field_animation_fname, 64, ImGuiInputTextFlags_CharsNoBlank);
    ImGui::PopItemWidth();


    /*****************************
     **   SIMULATION EXECUTION  **
     *****************************/

    // If the setup succeeds, this flag indicates that 
    // the simulation can be run
    static bool grid_ready = false;
    // Setup the grid according to the user input
    if (ImGui::Button("Setup Simulation")){
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
      for (size_t i = 0; i < g1->material_list.size(); i++){
        if (g1->material_list[i].x2 >= selected_Nx1D){
          g1->material_list[i].x2 = selected_Nx1D - 1;
        }
      }
      grid_ready = true;
    }
    if (grid_ready){
      ImGui::SameLine();
      if (ImGui::Button("Run Simulation")){
        g1->run_simulation();
      }
    }

    // The grid can be reset after it is created or the simulation run
    if (grid_ready){
      if (ImGui::Button("Reset grid")){
        delete g1;
        g1 = new Grid1D();
        grid_ready = false;
      }
    }

    ImGui::EndGroup();
    ImGui::EndGroup();

    /*****************************
     **   Material List Review  **
     *****************************/

    static bool show_mats_in_grid = false;
    if (ImGui::Button("Show materials")){
      show_mats_in_grid = show_mats_in_grid ? false : true;
    }

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

  return 0;
}

static void glfw_error_callback(int error, const char* description)
{
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}
