
expected_ants= c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP")

selected_colonies<-list(
  "B16-1", "B16-2", "B16-3", "B16-4","B16-5","B16-6","B16-7", "B16-8",
  "b16-1", "b16-2", "b16-3", "b16-4", "b16-5", "b16-6", "b16-7","b16-8",
  "bB16-1","bB16-2", "bB16-3", "bB16-4", "bB16-5","bB16-6","bB16-7","bB16-8",
  "bA16-1","bA16-2", "bA16-3", "bA16-4", "bA16-5","bA16-6","bA16-7","bA16-8"
)
#adjmatrix_list<-get_colony_files(source_folder_adjmatrix, selected_colonies)

present_ants_list<-list(  "B16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "B16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"), 
                          "B16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "B16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "B16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "B16-6" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "B16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "B16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          
                          "b16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "b16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "b16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "b16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "b16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "b16-6" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "b16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "b16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          
                          "bB16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bB16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bB16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bB16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bB16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bB16-6" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "PB", "PG", "PO", "PP"), #15 OP IS MISSING
                          "bB16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bB16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          
                          "bA16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bA16-2" = c("BB", "BG", "BO", "BP", "GB", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bA16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bA16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bA16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PO", "PP"),
                          "bA16-6" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bA16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                          "bA16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP")
)


framerate<-list(  "B16-1" =10,
                  "B16-2" = 10,
                  "B16-3" = 10,
                  "B16-4"= 10,
                  "B16-5"= 10,
                  "B16-6" = 10,
                  "B16-7" = 10,
                  "B16-8" = 10,
                  "b16-1"= 10,
                  "b16-2"= 10,
                  "b16-3" = 10,
                  "b16-4"= 10,
                  "b16-5"=10,
                  "b16-6" = 10,
                  "b16-7" = 10,
                  "b16-8"=10,
                  "bB16-1" =10,
                  "bB16-2"=10,
                  "bB16-3"=10,
                  "bB16-4"=10,
                  "bB16-5"=10,
                  "bB16-6"=10,
                  "bB16-7"=10,
                  "bB16-8"=10,
                  "bA16-1"=10,
                  "bA16-2"=10,
                  "bA16-3"=10,
                  "bA16-4"=10,
                  "bA16-5"=10,
                  "bA16-6"=10,
                  "bA16-7"=10,
                  "bA16-8"=10)


treatment_match<-list( "B16-1" = 1, 
                       "B16-2" = 1,
                       "B16-3" = 1,
                       "B16-4" = 1,
                       "B16-5" = 1,
                       "B16-6" = 1,
                       "B16-7" = 1,
                       "B16-8" = 1,
                       
                       "b16-1" = 2,
                       "b16-2" = 2,
                       "b16-3" = 2,
                       "b16-4" = 2,
                       "b16-5" = 2,
                       "b16-6" = 2,
                       "b16-7" = 2,
                       "b16-8" = 2,
                       
                       "bB16-1" = 3,
                       "bB16-2" = 3,
                       "bB16-3" = 3,
                       "bB16-4" = 3,
                       "bB16-5" = 3,
                       "bB16-6" = 3,
                       "bB16-7" = 3,
                       "bB16-8" = 3,
                       
                       "bA16-1" = 4,
                       "bA16-2" = 4,
                       "bA16-3" = 4,
                       "bA16-4" = 4,
                       "bA16-5" = 4,
                       "bA16-6" = 4,
                       "bA16-7" = 4,
                       "bA16-8" = 4
)
parameters <-c("total_number_interactions", "total_sum_of_interactions", "strength_mean", "strength_sd", "interaction_length_mean", "interaction_length_sd", "waiting_time_mean", "waiting_time_sd", "burstiness", "density_aggr", "global_eff", "assortativity", "mean_distance", "clustering_global")
