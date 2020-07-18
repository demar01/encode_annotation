final_annotated<-read_tsv(here("data/encode/k562/deseq_k562_annotated.txt"))

set.seed(123)
K562_split <- initial_split(final_annotated, strata = which_localization)
K562_train <- training(K562_split)
K562_test <- testing(K562_split)

K562_recipe <- recipe(which_localization ~ ., data = K562_train) %>%
  step_novel(all_predictors(), -all_numeric()) %>% 
  step_unknown(all_predictors(), -all_numeric()) %>% 
  step_medianimpute(all_numeric()) %>% 
  update_role(ensembl_gene_id_version, new_role = "id") %>% 
  step_dummy(all_predictors(), -all_numeric())


K562_wf <- workflow() %>%
  add_recipe(K562_recipe)


tree_spec <- bag_tree() %>% set_engine("rpart", times = 25) %>%
  set_mode("regression")

mars_spec <- bag_mars() %>%set_engine("earth", times = 25) %>%
  set_mode("regression")

tree_rs <- K562_wf %>%
  add_model(tree_spec) %>%
  fit(K562_train)

mars_rs <- K562_wf %>% add_model(mars_spec) %>%
  fit(K562_train)


test_rs <- K562_test %>%
  bind_cols(predict(tree_rs, K562_test)) %>%
  rename(.pred_tree = .pred) %>%
  bind_cols(predict(mars_rs, astro_test)) %>%
  rename(.pred_mars = .pred)

test_rs