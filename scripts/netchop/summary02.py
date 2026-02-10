import pandas as pd
import argparse
import os
import sys

def process_annotsv(file_path):
    """处理新的 AnnotSV 数据"""
    print(f"1. 读取 AnnotSV 数据: {file_path}")
    df = pd.read_csv(file_path)
    
    # 添加 Source_Tool 标签
    df['Source_Tool'] = 'annotsv_pep'
    
    # 确保 Mutation_Type 列存在 (如果还没有改名，尝试自动改)
    # 假设用户已经改好了，这里只做检查
    if 'Mutation_Type' not in df.columns:
        # 容错：如果还在 Source_Tool 列里但还没改名，这里就不自动处理了，
        # 因为脚本目的是合并。假设输入已经是清洗好的。
        # 如果真的缺，就填 Unknown
        print("   [提示] 未找到 Mutation_Type 列，填充为 'Unknown'")
        df['Mutation_Type'] = 'Unknown'

    return df

def process_neosv(file_path):
    """处理旧的 NeoSV 数据"""
    print(f"2. 读取并筛选 NeoSV 数据: {file_path}")
    df = pd.read_csv(file_path)
    
    # 筛选：只保留 Source_Tool 为 'neosv' 的行
    original_count = len(df)
    df = df[df['Source_Tool'] == 'neosv'].copy()
    print(f"   筛选: {original_count} -> {len(df)} 行 (只保留 neosv)")
    
    # 补齐缺失的列 (对齐 AnnotSV 的结构)
    df['Mutation_Type'] = 'NA' # NeoSV 结果里可能没这个信息
    df['sv_type'] = 'NA'       # NeoSV 结果里可能没这个信息
    
    return df

def main():
    parser = argparse.ArgumentParser(description="Merge AnnotSV (New) and NeoSV (Old) results.")
    parser.add_argument('--annot', required=True, help="新的 AnnotSV csv (summary.csv)")
    parser.add_argument('--neo', required=True, help="旧的 NeoSV/NetChop csv")
    parser.add_argument('--output', default="peptide_summary_final.csv", help="最终输出文件名")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.annot) or not os.path.exists(args.neo):
        print("错误：输入文件不存在。")
        sys.exit(1)
        
    # 处理
    df_1 = process_annotsv(args.annot)
    df_2 = process_neosv(args.neo)
    
    # 合并
    print("3. 执行合并...")
    final_df = pd.concat([df_1, df_2], ignore_index=True)
    
    # 整理列顺序 (美观)
    target_cols = [
        'Ident', 'peptide_sequence', 'peptide_length', 'sample', 
        'Gene', 'Source_Tool', 'Mutation_Type', 'sv_type', 
        'N_score', 'C_score', 'max_internal_score', 'internal_cleavage_count'
    ]
    
    # 确保只包含存在的列
    final_cols = [c for c in target_cols if c in final_df.columns]
    # 把其他可能存在的列放到后面
    remaining = [c for c in final_df.columns if c not in final_cols]
    final_df = final_df[final_cols + remaining]
    
    # 输出
    final_df.to_csv(args.output, index=False)
    print("-" * 30)
    print(f"✅ 合并完成！")
    print(f"输出文件: {args.output}")
    print(f"总行数: {len(final_df)}")
    print(f"来源分布:\n{final_df['Source_Tool'].value_counts()}")

if __name__ == "__main__":
    main()
    """
    输入画图函数之前，合并两个工具来源
    python summary02.py \
    --annot peptide_summary.csv \
    --neo /sibpt/anwy/lungfile/filter/netchop/peptide_summary.csv \
    --output peptide_summary_final.csv
        
    """