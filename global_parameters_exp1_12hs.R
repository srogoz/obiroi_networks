#################
#GLOBAL PARAMETERS EXP1_12hs

expected_ants= c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP")

selected_colonies<-list(
  "a16-1", "a16-2", "a16-3", "a16-4","a16-5","a16-6","a16-7", "a16-8",
  "b16-1", "b16-2", "b16-3", "b16-4", "b16-5", "b16-6", "b16-7","b16-8",
  "ba16-1","ba16-2", "ba16-3", "ba16-4", "ba16-5",  "ba16-6","ba16-7","ba16-8"
)

selected_colonies_dot<-list(
  "a16.1", "a16.2", "a16.3", "a16.4","a16.5","a16.6","a16.7", "a16.8",
  "b16.1", "b16.2", "b16.3", "b16.4", "b16.5", "b16.6", "b16.7","b16.8",
  "ba16.1","ba16.2", "ba16.3", "ba16.4", "ba16.5",  "ba16.6","ba16.7","ba16.8"
)

#EXP1 12hs
present_ants_list<-list( "a16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-6" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO"),
                         "b16-6" = c("BB", "BG", "BO", "BP", "GB", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OO", "PB", "PG", "PO", "PP"),
                         "ba16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-3" = c("BB", "BG", "BO", "BP", "GB", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-6" = c("BB", "BG", "BO", "BP", "GB", "GO", "OB", "OG", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-7" = c("BB", "BG", "GB", "GG", "GO", "GP", "OB", "OG", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP")
)

framerate<-list("a16-1" = 10,
                "a16-2" = 10,
                "a16-3" = 10,
                "a16-4" = 10,
                "a16-5" = 10,
                "a16-6" = 10,
                "b16-1" = 10,
                "b16-2" = 10,
                "b16-3" = 10,
                "b16-4" = 10,
                "b16-5" = 10,
                "b16-6" = 10,
                "ba16-1"= 10,
                "ba16-2"= 10,
                "ba16-3"= 10,
                "ba16-4"= 10, 
                "ba16-5"= 10,
                "ba16-6"= 10,
                "a16-7" = 5,
                "a16-8" = 5,
                "b16-7" = 5,
                "b16-8" = 5,
                "ba16-7"= 5,
                "ba16-8"= 5)


treatment_match<-list("a16-1" = "1",
                      "a16-2" = "1",
                      "a16-3" = "1",
                      "a16-4" = "1",
                      "a16-5" = "1",
                      "a16-6" = "1",
                      "a16-7" = "1",
                      "a16-8" = "1",
                      "b16-1" = "2", 
                      "b16-2" = "2",
                      "b16-3" = "2", 
                      "b16-4" = "2", 
                      "b16-5" = "2", 
                      "b16-6" = "2",
                      "b16-7" = "2",
                      "b16-8" = "2",
                      "ba16-1" = "3", 
                      "ba16-2" = "3",
                      "ba16-3" = "3",
                      "ba16-4" = "3",
                      "ba16-5" = "3",
                      "ba16-6" = "3",
                      "ba16-7" = "3",
                      "ba16-8" = "3"
)
genotype_match<-list("BB" = "b",
                     "BG" = "b",
                     "BO" = "b",
                     "BP" = "b",
                     "GB" = "b",
                     "GG" = "b", 
                     "GO" = "b",
                     "GP" = "b", 
                     "OB" = "a",
                     "OG" = "a",
                     "OO" = "a",
                     "OP" = "a",
                     "PB" = "a",
                     "PG" = "a",
                     "PO" = "a",
                     "PP" = "a")

exposed_ants <-c(
  "a16-1" = "OP",
  "a16-2" = "BP",
  "a16-3" = "OO",
  "a16-4" = "BO",
  "a16-5" = "OG",
  "a16-6" = "BG",
  "a16-7"= "OB",
  "a16-8"= "BB",
  "b16-1" = "PP",
  "b16-2" = "GP",
  "b16-3"= "PO",
  "b16-4" = "GO",
  "b16-5" = "PG",
  "b16-6" = "GG",
  "b16-7"= "PB",
  "b16-8"= "GB",
  "ba16-1"= "OP",
  "ba16-2"="PP",
  "ba16-3"="OG",
  "ba16-4"="GO", 
  "ba16-5"="GG",
  "ba16-6"="PB",
  "ba16-7"="GP",
  "ba16-8"="BB"
)


exposed_ants_string <-c(
  "OP_a16.1",
  "BP_a16.2",
  "OO_a16.3",
  "BO_a16.4",
  "OG_a16.5",
  "BG_a16.6",
  "OB_a16.7",
  "BB_a16.8",
  "PP_b16.1",
  "GP_b16.2",
  "PO_b16.3",
  "GO_b16.4",
  "PG_b16.5",
  "GG_b16.6",
  "PB_b16.7",
  "GB_b16.8",
  "OP_ba16.1",
  "PP_ba16.2",
  "OG_ba16.3",
  "GO_ba16.4", 
  "GG_ba16.5",
  "PB_ba16.6",
  "GP_ba16.7",
  "BB_ba16.8"
)

exposed_ants_plot <-c(
  "a16.1" ="OP",
  "a16.2" = "BP",
  "a16.3" = "OO",
  "a16.4" = "BO",
  "a16.5" = "OG",
  "a16.6" = "BG",
  "a16.7" = "OB",
  "a16.8" = "BB",
  "b16.1" = "PP",
  "b16.2" = "GP",
  "b16.3" = "PO",
  "b16.4" = "GO",
  "b16.5" = "PG",
  "b16.6" = "GG",
  "b16.7" = "PB",
  "b16.8" = "GB",
  "ba16.1" = "OP",
  "ba16.2" = "PP",
  "ba16.3" = "OG",
  "ba16.4" = "GO", 
  "ba16.5" = "GG",
  "ba16.6" = "PB",
  "ba16.7" = "GP",
  "ba16.8" = "BB"
)

parameters <-c("total_number_interactions", "total_sum_of_interactions", "strength_mean", "strength_sd", "interaction_length_mean", "interaction_length_sd", "waiting_time_mean", "waiting_time_sd", "burstiness", "density_aggr", "global_eff", "assortativity", "mean_distance", "clustering_global")

#plot colors
treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

anttypes_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536"
  
)

infection_status<-c("before" = "gray75",
                    "after" = "seagreen")

colony_colors <-c(
  "a16-1" = "blue",
  "a16-2" = "royalblue",
  "a16-3" = "slateblue2",
  "a16-4" = "navyblue",
  "a16-5" = "lightblue2",
  "a16-6" = "deepskyblue1",
  "a16-7"="skyblue2",
  "a16-8"="cyan",
  
  "b16-1" = "firebrick2",
  "b16-2" = "red4",
  "b16-3"= "tomato2",
  "b16-4" = "indianred",
  "b16-5" = "darkorange2",
  "b16-6" = "sienna1",
  "b16-7"="coral",
  "b16-8"="lightsalmon",
  
  "ba16-1"= "darkviolet",
  "ba16-2"="orchid2",
  "ba16-3"="darkmagenta",
  "ba16-4"="magenta1", 
  "ba16-5"="deeppink1",
  "ba16-6"="plum2",
  "ba16-7"="violet",
  "ba16-8"="hotpink2"
)

global_times<-c("0-1.25h", "1.25-2.5h","8.75h-10h", "10h-11.25h", "fungus exp")
global_times_colors_a<-c("0-1.25h" = "cyan", 
                         "1.25-2.5h" = "slateblue2",
                         "8.75h-10h" = "blue",
                         "10h-11.25h" = "dodgerblue1",
                         "fungus exp"= "springgreen1")


global_times_colors_b<-c("0-1.25h" = "firebrick2", 
                         "1.25-2.5h" = "red4",
                         "8.75h-10h" = "tomato2",
                         "10h-11.25h" = "darkorange",
                         "fungus exp"= "springgreen1")

global_times_colors_ba<-c("0-1.25h" = "darkviolet", 
                          "1.25-2.5h" = "magenta",
                          "8.75h-10h" = "plum2",
                          "10h-11.25h" = "orchid",
                          "fungus exp"= "springgreen1")



