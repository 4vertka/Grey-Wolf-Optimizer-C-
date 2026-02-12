
cmake_minimum_required(VERSION 3.16)                                    
project(project LANGUAGES CXX)                                          
                                                                        
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)                                   
                                                                        
add_executable(project main.cpp)                                        
                                                                        
target_link_libraries(project PUBLIC matplot)                           
