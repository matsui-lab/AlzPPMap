#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(igraph)
library(plotly)
library(shinybusy)
library(DT)
library(data.table)

# グラフをRDSファイルから読み込む
g <- readRDS("network_graph.rds")

# データの事前読み込みとキャッシュ
dist_data <- fread("eculidean_all.txt")
dist_data <- dist_data[, -1, with = FALSE] 

function(input, output, session) {
  
  gene_list <- reactive({
    rownames(as.matrix(get.adjacency(g)))
  })
  
  # 遺伝子選択のUIを生成
  output$gene_selector <- renderUI({
    selectizeInput("genes", "Select genes (up to 10):",
                   choices = NULL,
                   selected = NULL, 
                   multiple = TRUE,
                   options = list(maxItems = 10, 
                                  load = I("function(query, callback) { 
                                  Shiny.setInputValue('load_genes', query, {priority: 'event'}); 
                                }")))
  })
  
  # 初期選択とサーバーサイドの選択肢更新
  observe({
    updateSelectizeInput(session, "genes", choices = gene_list(), 
                         selected = c("HSPA5", "APP", "MAPT"), 
                         server = TRUE)
  })
  
  # スライダーの動的生成
  output$slider_inputs <- renderUI({
    lapply(input$genes, function(gene) {
      sliderInput(paste0("slider_", gene), paste("Select top % for", gene),
                  min = 1, max = 20, value = 5)
    })
  })
  
  # 分析実行のボタンアクション
  observeEvent(input$run_analysis, {
    req(input$genes)
    show_modal_spinner(spin = "circle", text = "Running analysis...")
    
    # 効率的なデータ処理: 選択された各遺伝子に対して、関連するデータを効率的に抽出
    get_top_percent_by_gene <- function(dist_data, gene_name, top_percent) {
      # data.tableを使用して効率的にデータをフィルタリング
      related_to_gene <- dist_data[Gene1 == gene_name | Gene2 == gene_name]
      setorder(related_to_gene, Distance)
      top_n <- floor(top_percent / 100 * nrow(related_to_gene))
      top_genes <- related_to_gene[1:top_n, ]
      return(unique(c(top_genes$Gene1, top_genes$Gene2)))
    }
    
    # 選択された全ての遺伝子に対して、上位の関連遺伝子を取得
    top_genes_list <- lapply(input$genes, function(gene) {
      get_top_percent_by_gene(dist_data, gene, input[[paste0("slider_", gene)]])
    })
    
    # 全ての遺伝子に共通する遺伝子を取得
    common_genes <- Reduce(intersect, top_genes_list)
    
    if (length(common_genes) == 0) {
      showModal(modalDialog(
        title = "No Common Genes Found",
        "No common genes were found for the selected criteria. Please adjust the selections and try again.",
        easyClose = TRUE,
        footer = NULL
      ))
      remove_modal_spinner()
      return(NULL)
    }
    
    # サブグラフの構築
    sub_vertices <- which(V(g)$name %in% common_genes)
    sub_g <- induced_subgraph(g, sub_vertices)
    
    deg <- degree(sub_g)
    isolated_nodes <- V(sub_g)[deg == 0]
    sub_g <- delete_vertices(sub_g, isolated_nodes)
    
    community <- cluster_leading_eigen(sub_g)
    V(sub_g)$community <- community$membership
    
    set.seed(1)
    layout <- layout_with_fr(sub_g)
    
    community_colors <- rainbow(length(unique(community$membership)))
    
    # プロットの生成
    output$network_plot <- renderPlotly({
      p <- plot_ly()
      
      for (i in unique(V(sub_g)$community)) {
        community_nodes <- which(V(sub_g)$community == i)
        p <- p %>%
          add_trace(
            type = 'scatter',
            mode = 'markers+text',
            x = layout[community_nodes, 1],
            y = layout[community_nodes, 2],
            marker = list(color = community_colors[i], size = 7,
                          line = list(color = 'black', width = 0)),
            text = V(sub_g)$name[community_nodes],
            textposition = 'top right',
            hoverinfo = 'text+name',
            name = paste("Community", i)
          )
      }
      
      edges <- get.edgelist(sub_g)
      if (nrow(edges) > 0) {
        edge_shapes <- list()
        for (i in 1:nrow(edges)) {
          edge_shapes[[i]] <- list(
            type = "line",
            line = list(color = "rgba(0, 0, 0, 0.2)", width = 0.5),
            x0 = layout[which(V(sub_g)$name == edges[i,1]), 1],
            y0 = layout[which(V(sub_g)$name == edges[i,1]), 2],
            x1 = layout[which(V(sub_g)$name == edges[i,2]), 1],
            y1 = layout[which(V(sub_g)$name == edges[i,2]), 2]
          )
        }
        p <- p %>%
          layout(
            title = "Subnetwork",
            showlegend = TRUE,
            shapes = edge_shapes,
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            margin = list(l = 0, r = 0, b = 0, t = 0),
            legend = list(y = 0.8)
          )
      }
      
      p
    })
    
    community_df <- data.frame(
      Gene = V(sub_g)$name,
      Community = V(sub_g)$community,
      Degree = degree(sub_g)
    )
    
    # コミュニティデータの表示
    output$community_table <- renderDT({
      datatable(community_df, rownames = FALSE)
    })
    
    # データフレームをCSVとしてダウンロード
    output$download_table <- downloadHandler(
      filename = function() { "community_data.csv" },
      content = function(file) {
        write.csv(community_df, file, row.names = FALSE)
      }
    )
    
    remove_modal_spinner()
  })
  
  output$graph_download <- downloadHandler(
    filename = function() { "network_graph.rds" },
    content = function(file) {
      saveRDS(g, file)
    }
  )
}

