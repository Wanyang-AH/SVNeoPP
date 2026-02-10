import pandas as pd
import os
import argparse
import sys

def get_sample_name(filename):
    """
    从文件名提取样本名。
    假设文件名格式为: l056_AnnotSV.peptides....
    策略: 取第一个下划线前的部分，或者根据具体命名习惯调整
    """
    base = os.path.basename(filename)
    # 尝试提取 l0xx
    if "_" in base:
        return base.split("_")[0]
    return base.split(".")[0]

def main():
    parser = argparse.ArgumentParser(description="Combine sample CSVs and add 'sample' column.")
    parser.add_argument('-i', '--input', required=True, nargs='+', help="Input CSV files (e.g., *.csv)")
    parser.add_argument('-o', '--output', required=True, help="Final output CSV filename")
    
    args = parser.parse_args()
    
    all_dfs = []
    
    print(f"--- 开始合并 {len(args.input)} 个文件 ---")
    
    for csv_file in args.input:
        if not os.path.exists(csv_file):
            print(f"[跳过] 文件不存在: {csv_file}")
            continue
            
        try:
            df = pd.read_csv(csv_file)
            
            # 1. 提取样本名
            sample_name = get_sample_name(csv_file)
            
            # 2. 添加 sample 列 (插入到 Gene 前面比较好看)
            # 如果已经存在但为空，则覆盖；如果不存在，则创建
            df['sample'] = sample_name
            
            print(f"读取: {os.path.basename(csv_file)} -> 样本名: {sample_name} (行数: {len(df)})")
            all_dfs.append(df)
            
        except Exception as e:
            print(f"[错误] 读取 {csv_file} 失败: {e}")

    if not all_dfs:
        print("没有数据被合并。")
        sys.exit(1)

    # 3. 合并
    final_df = pd.concat(all_dfs, ignore_index=True)
    
    # 4. 调整列顺序 (使其与旧脚本兼容，同时保留 sv_type)
    # 理想顺序: Ident, peptide_sequence, peptide_length, sample, Gene, Source_Tool, sv_type, Scores...
    target_cols = [
        'Ident', 'peptide_sequence', 'peptide_length', 'sample', 
        'Gene', 'Mutation_Type', 'sv_type', 
        'N_score', 'C_score', 'max_internal_score', 'internal_cleavage_count'
    ]
    
    # 只保留存在的列
    final_cols = [c for c in target_cols if c in final_df.columns]
    
    # 如果还有其他新列（比如之前可能存在的 extra info），也加上
    remaining_cols = [c for c in final_df.columns if c not in final_cols]
    final_cols.extend(remaining_cols)
    
    final_df = final_df[final_cols]
    
    # 5. 输出
    final_df.to_csv(args.output, index=False)
    
    print("-" * 30)
    print(f"✅ 合并完成！")
    print(f"输出文件: {args.output}")
    print(f"总行数: {len(final_df)}")
    print(f"包含样本: {final_df['sample'].unique()}")

if __name__ == "__main__":
    main()